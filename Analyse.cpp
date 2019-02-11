/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: <antoine.bois.1@ulaval.ca>
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <pthread.h>
#include <zlib.h>
#include <unordered_map>
#include <set>
#include <string>
#include <sys/stat.h>
#include "eval.h"

#if defined(_M_X64) || defined(__amd64__)
#define CONVERSION (unsigned long)
#else
#define CONVERSION (unsigned int)
#endif

#define BUFFER_SIZE 256  /* Grosseur des tampons */
#define BOOTSTRAP 10000  /* Nombre d'échantillons bootstrap */
#define SLEEP_TIME 5     /* Taux de rafraichissement du thread affichant la
			  * progression (en secondes) */
#define get_CI(std, nb_iters) (1.96 * std) / sqrt (nb_iters) /* Intervalle
							      * de confiance
							      */
using namespace std;

/*
 * Types des variables.
 */
enum {BOOLEAN, ACCUMUL, DISCRETE, CUSTOM_BOOLEAN, LOC_CALC, GLOB_CALC};

/*
 * Types des opérations (calculs).
 */
enum {AND, OR, EQ, NE, GT, GE, LT, LE};

/*
 * Sauvegarder la valeur d'une variable pour l'individu en cours.
 */
union last_value
{
  char      *string_value;
  double    num_value;
};

/*
 * Pour pouvoir afficher la progression.
 */
struct pBar_args
{
  FILE      *progress_file;
  int       *progress;
  int       total_iters;
  int       threads_count;
};

/*
 * Les structs qui suivent auraient pu grandement être simplifiés par
 * l'usage de variables globales. À modifier?
 *
 * On recopie toutes les infos à chaque fois au lieu de simplement
 * passer un pointeur vers le struct de configuration afin d'éviter un
 * double déférencement de pointeurs.
 */

/*
 * Un struct pour contenir les arguments passés aux threads (obligatoire).
 */
struct thread_args
{
  /* Expressions booléennes (custom bool) */
  char         ***data_comp_list;
  int          **comp_op_list;
  int          **bool_op_list;

  int          **c_bool_vars_ranks;
  int          *c_bool_vars_count;
  unsigned int *c_bool_results;
  int          c_bool_count;

  /* Calculs locaux et conditionnels */
  char         **loc_list;
  int          **loc_vars_ranks;
  int          *loc_vars_count;
  double       *loc_results;
  int          loc_count; /* loc = loc sans cond. */

  int          *cond_vars_rank;
  int          total_loc_count; /* total loc = loc + cond. */

  /* Variables booléennes */
  unsigned int *bool_results;
  int          bool_vars_count;

  /* Variables accumulatrices */
  double       *acc_results;
  int          acc_vars_count;

  /* Variables discrètes */
  unordered_map<string, unsigned int> *discrete_results;
  int          *vars_types;
  int          discrete_vars_count;

  /* Info sur le scénario à analyser */
  char         *name;
  char         *path;
  int          num_scen;
  int          lower_lim;
  int          upper_lim;

  /* Info générale */
  int          iters_count;
  int          vars_count;
  int          pop;
  int          *progress;
  int          progress_id;
};

/*
 * Un struct pour contenir la configuration désirée par l'utilisateur.
 */
struct conf_args
{
  /* Expressions booléennes */
  char        ***data_comp_list;
  char        **bool_labels;
  int         **bool_vars_ranks;
  int         **comp_op_list;
  int         **bool_op_list;
  int         *bool_vars_count;
  int         bool_count;

  /* Calculs globaux */
  char        **calcs_list;
  char        **calcs_labels;
  int         **calcs_relative_ranks;
  int         **calcs_vars_types;
  int         calcs_count;

  /* Calculs locaux et conditionnels */
  char        **loc_list;
  char        **loc_labels;
  int         **loc_vars_ranks;
  int         *loc_vars_count;
  int         *cond_vars_rank;
  int         loc_count;
  int         total_loc_count;

  /* ICER (Incremental cost-effectiveness ratio) */
  int         *ICR_vars_ranks;
  int         *ICR_vars_inv;
  int         *ICR_vars_types;
  int         ICR_vars_count;
  int         ICR_cmp_rank;
  int         ICR_cmp_type;

  /* Variables à ne pas afficher */
  int         *no_show;

  int         discrete_vars_count;
};

/*
 * Un struct pour contenir les informations nécessaires à l'affichage
 * des résultats.
 */
struct print_func_args
{
  /* Variables standards */
  char         **vars_list;
  int          *vars_types;

  unsigned int *bool_results;
  double       *acc_results;
  unordered_map<string, unsigned int> *discrete_results;
  int          acc_vars_count;
  int          discrete_vars_count;
  int          bool_vars_count;

  /* Calculs globaux */
  char         **calcs_list;
  char         **calcs_labels;
  int          **calcs_relative_ranks;
  int          **calcs_vars_types;
  int          calcs_count;

  /* Calculs locaux et conditionnels */
  char         **loc_list;
  char         **loc_labels;
  double       *loc_results;
  int          *cond_vars_rank;
  int          loc_count;
  int          total_loc_count;

  /* Expressions booléennes */
  char         **c_bool_labels;
  unsigned int *c_bool_results;
  int          c_bool_count;

  /* Ce que l'on ne veut pas afficher */
  int          *no_show;

  /* ICER */
  double       *cmp_means;
  double       *cmp_stds;
  double       *ICR_vars_means;
  double       *ICR_vars_stds;
  int          *ICR_vars_ranks;
  int          *ICR_vars_types;
  int          ICR_vars_count;
  int          cmp_rank;
  int          cmp_type;

  /* Info générale */
  char         *name;
  int          num_scen;
  int          iters_count;
  int          vars_count;
  int          pop;
};

void  parse_configuration    (FILE * pConf, conf_args * to_fill,
			      int *vars_types, char **vars_list,
			      int vars_count);

int   find_var_rank          (char *var, int vars_count, char **vars_list);

void  find_var_type_and_rank (char *var, int vars_count, char **vars_list,
			      int *vars_types, int calcs_count,
			      char **calcs_labels, int loc_count,
			      char **loc_labels, int bool_count,
			      char **bool_labels, int *type, int *rank);

void  check_delim_exist      (char delim);

void  *parse_csv             (void *ptr);

void  *display_progress      (void *ptr);

void  print_results          (print_func_args *args);

int   sort_func              (const void *elem1, const void *elem2);

int main (int argc, char **argv)
{
  /* On appelle le script de lancement qui réappelle le programme avec
   * le bon nombre d'arguments */
  if (argc < 5)
    {
      char cmd [1024] = "lancer_analyse.sh";
      for (int i = 1; i < argc; i++)
	{
	  strcat (cmd, " ");
	  strcat (cmd, argv[i]);
	}

      exit (system (cmd));
    }

  /* Erreurs possibles.. */
  else if (argc == 5)
    {
      puts("Aucun scénario à analyser.");
      exit(1);
    }

  int iters_count = atoi(argv[3]);
  if (!iters_count)
    {
      puts("Aucune simulation pour le(s) scénario(s).");
      exit(1);
    }

  /* Recevoir les scénarios. */
  int scenarios_count = argc - 5;

  /* Shallow copy des arguments en ligne de cmd */
  char **scenarios_list = (char**) malloc (sizeof(char*) * scenarios_count);
  for (int i = 0; i < scenarios_count; ++i)
    {
      scenarios_list[i] = argv[5+i];
    }

  /* Ouvrir un fichier Output.gz pour extraire certaines infos. */
  char file_path  [BUFFER_SIZE];
  strcpy (file_path, argv[1]);
  strcat (file_path, "/Results/");
  strcat (file_path, scenarios_list[0]);
  strcat (file_path, "/0_Output.gz");

  gzFile pFile = gzopen (file_path, "rb");
  if (pFile == Z_NULL)
    {
      printf("Mauvais chemin ou fichier inexistant: %s\n", file_path);
      exit(1);
    }

  /* Pour comparaison plus tard avec fichier binaire */
  struct stat modif_time_buff;

  stat (file_path, &modif_time_buff);
  time_t last_modif_Gz = modif_time_buff.st_mtime;

  char line [BUFFER_SIZE];

  gzgets (pFile, line, BUFFER_SIZE);
  gzgets (pFile, line, BUFFER_SIZE);

  int  vars_count = 0;
  char *pch       = line;

  /* Compter le nombre de variables */
  while (*pch != '\n')
    {
      if (*pch == ',')
	vars_count++;
      pch++;
    }

  if (!vars_count)
    {
      printf("Aucune variable détectée dans le fichier: %s\n", file_path);
      exit(1);
    }

  int *vars_types = (int*) malloc (sizeof(int) * vars_count);
  int bool_vars_count = 0;

  /* Les types sont trouvés automatiquement. Le type accumulateur
   * peut être changé en type discret par le fichier de conf.
   */
  strtok (line,",");
  for (int i = 0; i < vars_count; ++i)
    {
      pch = strtok (NULL, ",");
      if ( strncmp(pch, "true", 4) && strncmp(pch, "false", 5) )
	vars_types[i] = ACCUMUL ;
      else
	{
	  vars_types[i] = BOOLEAN ;
	  bool_vars_count++;
	}
    }

  gzclose (pFile);

  /* Ouvrir un fichier Summary.gz pour extraire certaines infos */
  strcpy(file_path+strlen(argv[1])+9+strlen(scenarios_list[0]),
	 "/0_Summary.gz");

  pFile = gzopen(file_path, "rb");
  if (pFile == Z_NULL)
    {
      printf("Mauvais chemin ou fichier inexistant: %s\n", file_path);
      exit(1);
    }

  char **vars_list = (char**) malloc (sizeof(char*) * vars_count);
  int  pop = 0 ;
  int  done_names_vars = 0 ;
  int  check_nb_vars ;

  /* La population: additionner les sous-populations. Les variables
   * doivent être les mêmes dans toutes les sous-populations: ceci est
   * est vérifié.  */
  while (gzgets(pFile, line, BUFFER_SIZE) != Z_NULL)
    {
      pch = strstr(line, "size=\"");
      if (pch != NULL)
	{
	  pop += atoi(strtok(pch+6, "\""));
	  check_nb_vars = 0;

	  /* Vrai seulement s'il s'agit d'une sous-population */
	  if (done_names_vars)
	    {
	      gzgets(pFile, line, BUFFER_SIZE);

	      while ( strstr (line, "</SubPopulation>") == NULL )
		{
		  strtok (line, "\"");
		  pch = strtok (NULL, "\"");

		  if ( strcmp (vars_list[check_nb_vars++], pch) )
		    {
		      printf ("La variable '%s' diffère de la population \
principale. Impossible d'effectuer l'analyse.", pch);
		      exit(1);
		    }

		  if (check_nb_vars > vars_count)
		    {
		      puts("Trop de variables dans une des \
sous-populations.");
		      exit(1);
		    }

		  gzgets(pFile, line, BUFFER_SIZE);
		}

	      if (check_nb_vars != vars_count)
		{
		  puts("Pas assez de variables dans une des \
sous-populations.");
		  exit(1);
		}
	    }
	  else
	    {
	      /* Population principale: on enregistre le nom des variables
	       * en mémoire. */
	      for (int i = 0; i < vars_count; i++)
		{
		  gzgets(pFile, line, BUFFER_SIZE);

		  strtok (line, "\"");
		  pch = strtok (NULL, "\"");
		  vars_list[i] = (char*) malloc (strlen(pch)+1);
		  strcpy (vars_list[i], pch);
		}
	      done_names_vars = 1;
	    }
	}
    }

  gzclose (pFile);

  if (!pop)
    {
      printf("Aucune population détectée dans le fichier: %s\n", file_path);
      exit(1);
    }

  /* Configuration de l'usager */
  FILE * pConf = fopen (argv[2], "r");
  conf_args current_conf;

  if (pConf != NULL)
    {
      parse_configuration (pConf, &current_conf, vars_types, vars_list,
			   vars_count);
      fclose (pConf);
    }
  else
    {
      /* Valeurs par défaut: ceci ne devrait jamais arriver car on crée
       * normalement un template vide si aucun fichier de configuration
       * n'est présent, qui pourra quand même être ouvert et analyser. */
      current_conf.discrete_vars_count = 0;
      current_conf.calcs_count         = 0;
      current_conf.loc_count           = 0;
      current_conf.total_loc_count     = 0;
      current_conf.bool_count          = 0;
      current_conf.ICR_vars_count      = 0;
      current_conf.no_show             = (int*) calloc (vars_count,
							sizeof(int) );
    }

  /* Pour éviter d'avoir à recalculer ces valeurs dans les ICR: tables
   * de correspondance selon le principe de mémoization */
  double *cmp_means = (double*) malloc (sizeof(double) * scenarios_count);
  double *cmp_stds  = (double*) malloc (sizeof(double) * scenarios_count);

  double *ICR_vars_means = (double*) malloc
    (sizeof(double) * scenarios_count * current_conf.ICR_vars_count);
  double *ICR_vars_stds  = (double*) malloc
    (sizeof(double) * scenarios_count * current_conf.ICR_vars_count);

  int acc_vars_count = vars_count - current_conf.discrete_vars_count
    - bool_vars_count;

  /* Alloué sur le monceau pour ne pas risquer de faire exploser la pile
   * dans les  grosses simulations (plusieurs scénarios, plusieurs
   * variables, etc.)
   */
  double *acc_results = (double*) malloc (sizeof(double) * scenarios_count
					  * acc_vars_count * iters_count);

  unsigned int *bool_results = (unsigned int*) malloc (sizeof(unsigned int)
     * bool_vars_count * scenarios_count * iters_count);

  unordered_map <string, unsigned int> *discrete_vars =
    new unordered_map <string, unsigned int>
    [scenarios_count * iters_count * current_conf.discrete_vars_count];

  double *loc_results = (double*) malloc (sizeof(double) * scenarios_count
     * current_conf.total_loc_count * iters_count);

  unsigned int *c_bool_results = (unsigned int*)
    malloc (sizeof(unsigned int) * current_conf.bool_count
	    * scenarios_count * iters_count);

  /* Pour le passage d'arguments à la fonction qui affiche les résultats */
  print_func_args *print_args = (print_func_args*) malloc
    (sizeof(print_func_args) * scenarios_count);

  print_args[0].vars_list            = vars_list;
  print_args[0].vars_types           = vars_types;
  print_args[0].name                 = scenarios_list[0];

  print_args[0].calcs_list           = current_conf.calcs_list;
  print_args[0].calcs_labels         = current_conf.calcs_labels;
  print_args[0].calcs_relative_ranks = current_conf.calcs_relative_ranks;
  print_args[0].calcs_vars_types     = current_conf.calcs_vars_types;

  print_args[0].bool_results         = bool_results;
  print_args[0].acc_results          = acc_results;
  print_args[0].discrete_results     = discrete_vars;

  print_args[0].cmp_means            = cmp_means;
  print_args[0].cmp_stds             = cmp_stds;
  print_args[0].ICR_vars_means       = ICR_vars_means;
  print_args[0].ICR_vars_stds        = ICR_vars_stds;
  print_args[0].ICR_vars_ranks       = current_conf.ICR_vars_ranks;
  print_args[0].ICR_vars_types       = current_conf.ICR_vars_types;
  print_args[0].cmp_rank             = current_conf.ICR_cmp_rank;
  print_args[0].cmp_type             = current_conf.ICR_cmp_type;
  print_args[0].calcs_count          = current_conf.calcs_count;

  print_args[0].num_scen             = 0;
  print_args[0].pop                  = pop;
  print_args[0].ICR_vars_count       = current_conf.ICR_vars_count;
  print_args[0].iters_count          = iters_count;
  print_args[0].vars_count           = vars_count;

  print_args[0].acc_vars_count       = acc_vars_count;
  print_args[0].discrete_vars_count  = current_conf.discrete_vars_count;
  print_args[0].bool_vars_count      = bool_vars_count;
  print_args[0].loc_count            = current_conf.loc_count;

  print_args[0].loc_list             = current_conf.loc_list;
  print_args[0].loc_labels           = current_conf.loc_labels;
  print_args[0].loc_results          = loc_results;

  print_args[0].c_bool_count         = current_conf.bool_count;
  print_args[0].c_bool_labels        = current_conf.bool_labels;
  print_args[0].c_bool_results       = c_bool_results;

  print_args[0].no_show              = current_conf.no_show;

  print_args[0].total_loc_count      = current_conf.total_loc_count;
  print_args[0].cond_vars_rank       = current_conf.cond_vars_rank;

  for (int i = 1; i < scenarios_count; i++)
    {
      print_args[i] = print_args[0]; /* deep copy */

      print_args[i].num_scen = i;
      print_args[i].name     = scenarios_list[i];
    }

  /* Vérifier si on a déjà fait le parsing avec les mêmes options.
   * Si oui, charger ces données en mémoire. Sinon, faire le parsing.
   */
  int *vars_types_check = (int*) malloc (sizeof(int) * vars_count);
  int same = 0;

  /* Le nom du fichier binaire: nom du fichier de configuration avec son
   * extension remplacée par '.aux' */
  char bin_path [BUFFER_SIZE];
  strcpy (bin_path, argv[1]);
  strcat (bin_path, "Analyse/");

  pch = argv[2];

  while (*pch)
    pch++;

  while (*pch != '/' && *pch != '.' && pch != argv[2])
    pch--;

  if (*pch == '.')
    *pch = NUL;

  while (*pch != '/' && pch != argv[2])
    pch--;

  strcat (bin_path, pch);
  strcat (bin_path, ".aux");

  FILE * pBin = fopen ( bin_path , "rb");

  /* Vérifie le nombre de lectures faites sur le fichier .aux */
  unsigned int read_count ;
  int keys_read_count = 0 ;
  int arg_counter     = 0 ;

  /* Vérifie qu'il n'y a pas de dépassements de tampon lors de la
   * vérification de calculs/expressions */
  int buff_overflow = 0;
  int finished      = 0; /* remplace des gotos (qui aurait pu être
			  * laissés et éviteraient des vérif. inutiles) */
  int i, v;

  if (pBin != NULL)
    {
      /* Fichier binaire existe: vérification que rien n'a changé. Voir
       * comment se fait la sérialization plus loins pour comprendre
       * comment le fichier binaire est décrypté.*/
      same = 1;
      read_count = fread (vars_types_check, sizeof(int), vars_count, pBin);

      int it_nb, loc_nb, c_bool_nb;
      read_count += fread (&it_nb, sizeof(int), 1 , pBin);
      read_count += fread (&loc_nb, sizeof(int), 1, pBin);
      read_count += fread (&c_bool_nb, sizeof(int), 1, pBin);

      if (loc_nb != current_conf.total_loc_count ||
	  c_bool_nb != current_conf.bool_count ||
          it_nb != iters_count)
	{
	  same = 0;
	  finished = 1;
	}

      int str_size;
      /* Calculs locaux */
      for (i = 0; i < current_conf.loc_count && !finished; i++)
	{
	  read_count += fread (&str_size, sizeof(int), 1, pBin);
	  if ( str_size > BUFFER_SIZE )
	    {
	      buff_overflow = finished = 1;
	      break;
	    }

	  read_count += fread (line, str_size, 1, pBin);

	  if ( strncmp (line, current_conf.loc_list[i], BUFFER_SIZE))
	    {
	      same = 0;
	      finished = 1;
	      break;
	    }
	}

      /* Calculs conditionnels */
      int cond_var_rank;
      for (i = current_conf.loc_count; i < current_conf.total_loc_count
	     && !finished;
	   i++)
	{
	  read_count += fread (&str_size, sizeof(int), 1, pBin);
	  if ( str_size > BUFFER_SIZE )
	    {
	      buff_overflow = finished = 1;
	      break;
	    }

	  read_count += fread (line, str_size, 1, pBin);
	  read_count += fread (&cond_var_rank, sizeof(int), 1, pBin);

	  if (strncmp(line, current_conf.loc_list[i], BUFFER_SIZE)
	      || cond_var_rank !=
	      current_conf.cond_vars_rank[i - current_conf.loc_count])
	    {
	      same = 0;
	      finished = 1;
	      break;
	    }
	}

      int arg_nb;
      int oper;
      /* Expressions booléennes  */
      for (i = 0; i < current_conf.bool_count && !finished; i++)
	{
	  read_count += fread (&arg_nb, sizeof(int), 1, pBin);
	  arg_counter += arg_nb;

	  if ( arg_nb != current_conf.bool_vars_count[i] )
	    {
	      same = 0;
	      finished = 1;
	      break;
	    }

	  for (v = 0; v < arg_nb && !finished; v++)
	  {
	    read_count += fread (&str_size, sizeof(int), 1, pBin);
	    if ( str_size > BUFFER_SIZE )
	      {
		buff_overflow = finished = 1;
		break;
	      }

	    read_count += fread (line, str_size, 1, pBin);
	    read_count += fread (&oper, sizeof(int), 1, pBin);

	    if ( strncmp (line, current_conf.data_comp_list[i][v],
			  BUFFER_SIZE)
		 || oper != current_conf.comp_op_list[i][v] )
	      {
		same = 0;
		finished = 1;
		break;
	      }
	  }

	  for (v = 0; v < arg_nb - 1 && !finished; v++)
	    {
	      read_count += fread(&oper, sizeof(int), 1, pBin);

	      if (oper != current_conf.bool_op_list[i][v])
		{
		  same = 0;
		  finished = 1;
		  break;
		}
	    }
	}

      /* Variables discrètes (indirectement) */
      for (i = 0; i < vars_count && !finished; i++)
	{
	  if (vars_types_check[i] != vars_types[i])
	    {
	      same = 0;
	      finished = 1;
	      break;
	    }
	}

      /* Si l'utilisateur s'amuse à copier le dossier 'Analyse' d'un projet
       * à un autre.. */
      stat (bin_path, &modif_time_buff);

      if (last_modif_Gz > modif_time_buff.st_mtime)
	same = 0;
    }

  free (vars_types_check);

  if (buff_overflow)
    {
      printf("Fichier binaire '%s' corrompu. Il est conseillé de le \
supprimer et de relancer ce programme.\n", bin_path);
      exit(1);
    }

  else if (same)
    {
      /* Les options n'ont pas changées et le parsing a déjà été fait:
       * charger les données en mémoire */
      read_count += fread (acc_results, sizeof(double), scenarios_count
			   * iters_count * acc_vars_count, pBin);

      read_count += fread (bool_results, sizeof(unsigned int),
			   scenarios_count * iters_count * bool_vars_count,
			   pBin);

      read_count += fread (loc_results, sizeof(double), scenarios_count
			   * iters_count * current_conf.total_loc_count,
			   pBin);

      read_count += fread (c_bool_results, sizeof(unsigned int),
			   scenarios_count * iters_count
			   * current_conf.bool_count, pBin);

      int keys_count;
      string current_key;
      int current_key_size;
      int value;

      for (i = 0; i < scenarios_count * iters_count
	     * current_conf.discrete_vars_count; ++i)
	{
	  read_count += fread (&keys_count, sizeof(int), 1, pBin);
	  keys_read_count += keys_count;

	  for (v = 0; v < keys_count; ++v)
	    {
	      read_count += fread (&current_key_size, sizeof(int), 1, pBin);
	      if ( current_key_size > BUFFER_SIZE )
		{
		  printf("Fichier binaire '%s' corrompu. Il est conseillé \
de le supprimer et de relancer ce programme.\n", bin_path);
		  exit(1);
		}

	      read_count += fread (line, current_key_size, 1, pBin);
	      read_count += fread (&value, sizeof(unsigned int), 1, pBin);
	      discrete_vars[i][string(line)] = value;
	    }
	}

      fclose(pBin);

      /* Sorte de checksum */
      if (read_count !=
	  3 + (current_conf.loc_count * 2) + arg_counter + vars_count
	  + ( (current_conf.total_loc_count - current_conf.loc_count +
	     keys_read_count + arg_counter) * 3 )
	  + ( scenarios_count * iters_count *
	      (acc_vars_count + bool_vars_count +
	       current_conf.total_loc_count + current_conf.bool_count +
	       current_conf.discrete_vars_count) ))
	{
	  puts("Certaines données n'ont pu être correctement récupérées. \
Cela peut être dû à un fichier '.aux' corrompu, des droits de lecture \
manquants, ou un système de fichier instable.");
	  exit(1);
	}

      /* Afficher les données qui furent chargées */
      for (i = 0; i < scenarios_count; ++i)
	print_results (print_args + i);
    }

  else
    {
      /* Parsing doit se faire! */
      if (pBin != NULL)
	fclose(pBin);

      /* Fichier qui montre la progression en temps réel */
      char progress_file [BUFFER_SIZE];
      strcpy (progress_file, argv[1]);
      strcat (progress_file, "progression.txt");

      int max_threads = atoi(argv[4]);

      /* Nombre de threads par scénario */
      int threads_per_scen =  (iters_count >= ( max_threads /
						      scenarios_count ) ) ?
	(int)ceil((float)max_threads / scenarios_count) : 1;

      /* La progression de chaque thread */
      int *progress = (int*) calloc (threads_per_scen
				     * scenarios_count, sizeof(int));

      /* Contenir les threads eux-mêmes */
      pthread_t *threads_array = (pthread_t*) malloc
	(sizeof(pthread_t) * scenarios_count * threads_per_scen);

      /* Contenir les paramètres à passer aux threads */
      thread_args *list_args = (thread_args*) malloc
	(sizeof(thread_args) * scenarios_count * threads_per_scen);

      /* Les paramètres à passer aux threads, en commençant par le premier
       * thread. */
      strcpy (file_path, argv[1]);
      strcat (file_path, "/Results/");

      list_args[0].lower_lim           = 0;
      list_args[0].upper_lim           = iters_count/threads_per_scen;
      list_args[0].iters_count         = iters_count;

      list_args[0].acc_results         = acc_results;
      list_args[0].bool_results        = bool_results;
      list_args[0].discrete_results    = discrete_vars;
      list_args[0].c_bool_results      = c_bool_results;
      list_args[0].loc_results         = loc_results;

      list_args[0].vars_types          = vars_types;
      list_args[0].vars_count          = vars_count;
      list_args[0].discrete_vars_count = current_conf.discrete_vars_count;
      list_args[0].acc_vars_count      = acc_vars_count;
      list_args[0].bool_vars_count     = bool_vars_count;

      list_args[0].pop                 = pop;
      list_args[0].num_scen            = 0;
      list_args[0].name                = scenarios_list[0];
      list_args[0].path                = file_path;

      list_args[0].loc_count           = current_conf.loc_count;
      list_args[0].loc_vars_count      = current_conf.loc_vars_count;
      list_args[0].loc_vars_ranks      = current_conf.loc_vars_ranks;
      list_args[0].loc_list            = current_conf.loc_list;

      list_args[0].c_bool_count        = current_conf.bool_count;
      list_args[0].c_bool_vars_count   = current_conf.bool_vars_count;
      list_args[0].c_bool_vars_ranks   = current_conf.bool_vars_ranks;
      list_args[0].comp_op_list        = current_conf.comp_op_list;
      list_args[0].bool_op_list        = current_conf.bool_op_list;
      list_args[0].data_comp_list      = current_conf.data_comp_list;

      list_args[0].cond_vars_rank      = current_conf.cond_vars_rank;
      list_args[0].total_loc_count     = current_conf.total_loc_count;

      list_args[0].progress            = progress;
      list_args[0].progress_id         = 0;

      /* Premier scénario */
      pthread_create(threads_array, NULL, parse_csv, (void*) list_args);

      for (v = 1; v < threads_per_scen; ++v)
	{
	  list_args[v]             = list_args[v-1]; /* deep copy */
	  list_args[v].lower_lim   = list_args[v-1].upper_lim;
	  list_args[v].progress_id = v;

	  list_args[v].upper_lim = (v == threads_per_scen - 1) ?
	    iters_count : list_args[v-1].upper_lim +
	    (list_args[v-1].upper_lim - list_args[v-1].lower_lim);

	  pthread_create(threads_array+v, NULL, parse_csv,
			 (void*) (list_args+v));
	}

      int offset;
      /* Les autres scénarios */
      for (i = 1; i < scenarios_count; ++i)
	{
	  offset            = i * threads_per_scen;
	  list_args[offset] = list_args[0];

	  list_args[offset].name        = scenarios_list[i];
	  list_args[offset].num_scen    = i;
	  list_args[offset].progress_id = offset;

	  if (threads_per_scen > 1)
	    {
	      list_args[offset].lower_lim = 0;
	      list_args[offset].upper_lim = iters_count /
		threads_per_scen;
	    }

	  pthread_create(threads_array + offset, NULL, parse_csv,
			 (void*) (list_args + offset));

	  /* Division des threads pour un même scénario */
	  for (v = 1; v < threads_per_scen; ++v)
	    {
	      list_args[offset+v] = list_args[offset+v-1];
	      list_args[offset+v].progress_id++;

	      list_args[offset+v].lower_lim = list_args[offset+v].upper_lim;

	      list_args[offset+v].upper_lim = (v == threads_per_scen
					       - 1) ? iters_count :
		list_args[offset+v].upper_lim
		+ (list_args[offset+v-1].upper_lim
		   - list_args[offset+v-1].lower_lim);

	      pthread_create(threads_array+offset+v, NULL, parse_csv,
			     (void*) (list_args+offset+v));
	    }
	}

      /* On crée le thread qui servira à suivre la progression */
      FILE *p_bar = fopen(progress_file, "w");
      pthread_t progress_thread;
      pBar_args current_progress = { p_bar, progress,
				     scenarios_count * iters_count,
				     threads_per_scen
				     * scenarios_count };

      pthread_create (&progress_thread, NULL, display_progress,
		      (void*) (&current_progress));

      /* Tous les threads ont été créés. Il faut traiter leur mort. */
      void *status;
      int  *done_parsing = (int*) calloc (scenarios_count, sizeof(int));

      for (i = 0; i < scenarios_count * threads_per_scen; ++i)
	{
	  pthread_join (threads_array[i], &status);
	  done_parsing [CONVERSION status] ++ ;

	  /* Si vrai, alors impossible qu'un scénario ait eu le temps
	   * de finir: on saute la vérification */
	  if (i < (threads_per_scen - 1))
	    {
	      continue;
	    }
	  else
	    {
	      for (v = 0; v < scenarios_count; ++v)
		{
		  if (done_parsing [v] == threads_per_scen)
		    {
		      /* Affichage des résultats */
		      print_results (print_args + v);
		      done_parsing [v] = -1;
		    }
		}
	    }
	}
      free (threads_array);
      free (list_args);
      free (done_parsing);

      /* Le parsing est complété, alors le thread de progression devient
       * inutile: il faut l'interrompre. */
      pthread_cancel (progress_thread);
      pthread_join   (progress_thread, NULL);
      fclose (p_bar);
      remove (progress_file);
      free (current_progress.progress);

      /* Sauvegarder les résultats du parsing */
      pBin = fopen ( bin_path , "wb");
      if (pBin != NULL)
	{
	  /* Même ordre que la lecture du fichier binaire plus haut */
	  fwrite (vars_types, sizeof(int), vars_count, pBin);

	  fwrite (&iters_count, sizeof(int), 1, pBin);

          fwrite (&current_conf.total_loc_count, sizeof(int), 1, pBin);

	  fwrite (&current_conf.bool_count, sizeof(int), 1, pBin);

	  int  key_size;
	  /* Calculs locaux et condtionnels */
	  for (i = 0; i < current_conf.total_loc_count; ++i)
	    {
	      key_size = strlen(current_conf.loc_list[i]) + 1;
	      fwrite (&key_size, sizeof(int), 1, pBin);
	      fwrite (current_conf.loc_list[i], key_size, 1, pBin);

	      if (i >= current_conf.loc_count)
		fwrite (&current_conf.cond_vars_rank
			[i - current_conf.loc_count],
			sizeof(int), 1 , pBin);
	    }

	  /* Expressions booléennes */
	  for (i = 0; i < current_conf.bool_count; ++i)
	    {
	      fwrite (&current_conf.bool_vars_count[i], sizeof(int),
		      1, pBin);

	      for (v = 0; v < current_conf.bool_vars_count[i]; v++)
		{
		  key_size = strlen(current_conf.data_comp_list[i][v]) + 1;
		  fwrite (&key_size, sizeof(int), 1, pBin);
		  fwrite (current_conf.data_comp_list[i][v], key_size,
			  1, pBin);
		  fwrite (&current_conf.comp_op_list[i][v], sizeof(int),
			  1, pBin);
		}

	      for (v = 0; v < current_conf.bool_vars_count[i] - 1; v++)
		{
		  fwrite (&current_conf.bool_op_list[i][v], sizeof(int),
			  1, pBin);
		}
	    }

	  /* Les résultats eux-mêmes */
	  fwrite (acc_results, sizeof(double), scenarios_count
		  * iters_count * acc_vars_count, pBin);

	  fwrite (bool_results, sizeof(unsigned int), scenarios_count
		  * iters_count * bool_vars_count, pBin);

	  fwrite (loc_results, sizeof(double), scenarios_count
		  * iters_count * current_conf.total_loc_count, pBin);

	  fwrite (c_bool_results, sizeof(unsigned int), scenarios_count
		  * iters_count * current_conf.bool_count, pBin);

	  int  keys_count;
	  const char *c_key;

	  /* Les résultats des variables discrètes est un peux plus
	   * complexe: il s'agit de séraliser un std::map */
	  for (i = 0; i < scenarios_count * iters_count
		 * current_conf.discrete_vars_count; ++i)
	    {
	      keys_count = 0;

	      for (unordered_map<string, unsigned int>::iterator it =
		     discrete_vars[i].begin(); it!=discrete_vars[i].end();
		   ++it)
		{
		  keys_count++;
		}

	      fwrite (&keys_count, sizeof(int), 1, pBin);

	      for (unordered_map<string, unsigned int>::iterator it =
		     discrete_vars[i].begin(); it!=discrete_vars[i].end();
		   ++it)
		{
		  /* Strings C plus facile à sérialiser */
		  c_key = it->first.c_str();
		  key_size = strlen(c_key) + 1;
		  fwrite (&key_size, sizeof(int), 1, pBin);

		  fwrite (c_key, key_size, 1, pBin);
		  fwrite (&it->second, sizeof(unsigned int), 1, pBin);
		}
	    }
	  fclose (pBin);
	}
      else
	{
	  printf("Incapable d'écrire le fichier binaire: %s\n\n", bin_path);
	  /* N'interrompt pas le programme mais devrait peut-être */
	}
    }

  /* Si pas de variables d'ICR, fin du programme */
  if (! current_conf.ICR_vars_count)
    return 0;

  /* Version très beta des incertitudes sur les ICER: il est très
   * encourager de trouver une meilleure méthode.
   */

  puts("* Calculs d\'ICER (Incremental Cost-Effectiveness Ratio)");
  puts("( À noter: IC = intervalle de confiance à ~ 95% )");

  /* On veut établir le rang des scénarios selon le comparateur */
  double *sorted_cmp_means = (double*) malloc
    (sizeof(double) * scenarios_count);

  memcpy (sorted_cmp_means, cmp_means, sizeof(double) * scenarios_count);
  qsort  (sorted_cmp_means, scenarios_count, sizeof(double), sort_func);

  /* Les rangs sont trouvés et stockés ici */
  int *scenarios_new_ranks = (int*) malloc (sizeof(int) * scenarios_count);
  for (i = 0; i < scenarios_count; ++i)
    {
      for (v = 0; v < scenarios_count; ++v)
	{
	  if (sorted_cmp_means[i] == cmp_means[v])
	    {
	      scenarios_new_ranks[i] = v;
	      break;
	    }
	}
    }

  int previous, how_much_backward; /* scénario de comparaison & offset dû
				    * aux scénarios dominés */

  double c_cmp, inter_c_cmp, c_var, inter_c_var, delta_cmp, delta_var;
  double p_cmp, inter_p_cmp, p_var, inter_p_var, ICR;

  double denom; /* scindé la division en deux */
  double *ICR_bootstrap = (double*) malloc (sizeof(double) * BOOTSTRAP);

  char *to_display; /* variable en traitement */

  /* Le seed du générateur de nombres réels aléatoires */
  srand48(time(NULL));

  /* Commencer par formater un peu les tableaux */
  for (i = 0; i < current_conf.ICR_vars_count; ++i)
    {
      /* Comparateur */
      switch (current_conf.ICR_cmp_type)
	{
	case BOOLEAN: case ACCUMUL:
	  to_display = vars_list [current_conf.ICR_cmp_rank];
	  break;

	case CUSTOM_BOOLEAN:
	  to_display = current_conf.bool_labels[current_conf.ICR_cmp_rank];
	  break;

	case LOC_CALC:
	  to_display = current_conf.loc_labels[current_conf.ICR_cmp_rank];
	  break;

	case GLOB_CALC:
	  to_display = current_conf.calcs_labels[current_conf.ICR_cmp_rank];
	  break;
	}

      printf("\nOption,Total %s,IC(±),Delta %s,", to_display, to_display);

      /* Variable au dénominateur */
      switch (current_conf.ICR_vars_types[i])
	{
	case BOOLEAN: case ACCUMUL:
	  to_display = vars_list[current_conf.ICR_vars_ranks[i]];
	  break;

	case CUSTOM_BOOLEAN:
	  to_display = current_conf.bool_labels
	    [current_conf.ICR_vars_ranks[i]];
	  break;

	case LOC_CALC:
	  to_display = current_conf.loc_labels[current_conf.ICR_vars_ranks
					       [i]];
	  break;

	case GLOB_CALC:
	  to_display = current_conf.calcs_labels
	    [current_conf.ICR_vars_ranks[i]];
	  break;
	}

      printf("Total %s,IC(±),Delta %s,ICR,IC(-),IC(+),Status\n",
	     to_display, to_display ) ;

      printf("%s,%.8G,%.8G,N/A,%.8G,%.8G,N/A,Baseline,N/A,N/A, \n",
	     scenarios_list[*scenarios_new_ranks], *sorted_cmp_means,
	     get_CI (cmp_stds[*scenarios_new_ranks], iters_count),
	     ICR_vars_means[*scenarios_new_ranks*current_conf
			    .ICR_vars_count+i],
	     get_CI (ICR_vars_stds[*scenarios_new_ranks*current_conf
				   .ICR_vars_count+i], iters_count));

      how_much_backward = 0; /* si scénario dominé, alors comparaison
			      * du suivant est faite avec le précédent */
      for (v = 1; v < scenarios_count; ++v)
	{
	  previous = v - 1 - how_much_backward;

	  /* Noms de variables plus clairs: c = current, p = previous */
	  c_cmp       = sorted_cmp_means[v];
	  inter_c_cmp = get_CI (cmp_stds[scenarios_new_ranks[v]]
				, iters_count);

	  p_cmp       = sorted_cmp_means[previous];
	  inter_p_cmp = get_CI (cmp_stds[scenarios_new_ranks[previous]]
				, iters_count);

	  delta_cmp   = c_cmp - p_cmp;

	  c_var       = ICR_vars_means[scenarios_new_ranks[v]
				       * current_conf.ICR_vars_count + i];
	  inter_c_var = get_CI (
			 ICR_vars_stds[scenarios_new_ranks[v]
				       * current_conf.ICR_vars_count + i]
			 , iters_count);

	  p_var       = ICR_vars_means[scenarios_new_ranks[previous]
				       * current_conf.ICR_vars_count + i];
	  inter_p_var = get_CI (
			 ICR_vars_stds[scenarios_new_ranks[previous]
				       * current_conf.ICR_vars_count + i]
			 , iters_count);

	  delta_var   = current_conf.ICR_vars_inv[i] ? p_var - c_var :
	    c_var - p_var;

	  if (! delta_var) /* div par zéro */
	    {
	      printf("%s,%.8G,%.8G,%.8G,%.8G,%.8G,%.8G,inf,N/A,N/A,Dominé\n"
		     ,scenarios_list[scenarios_new_ranks[v]], c_cmp,
		     inter_c_cmp, delta_cmp, c_var, inter_c_var, delta_var);
	      ++how_much_backward;
	      continue;
	    }

	  ICR = delta_cmp / delta_var;

	  /* On sélectionne des échantillons ici pour le calcul
	   * d'intervalle de confiance de l'ICR */
	  for (int k = 0; k < BOOTSTRAP; ++k)
	    {
	      denom = current_conf.ICR_vars_inv[i] ?
		((drand48() * inter_p_var * 2) + (p_var - inter_p_var))
		- ((drand48() * inter_c_var * 2) + (c_var - inter_c_var)):
		((drand48() * inter_c_var * 2) + (c_var - inter_c_var))
		- ((drand48() * inter_p_var * 2) + (p_var - inter_p_var));

	      ICR_bootstrap[k] = (( (drand48() * inter_c_cmp * 2)
				   + (c_cmp - inter_c_cmp) )
		- ( (drand48() * inter_p_cmp * 2)
		    + (p_cmp - inter_p_cmp) )) / denom;
	    }

	  /* On ordonne les ICR générés aléatoirement par bootstrap */
	  qsort (ICR_bootstrap, BOOTSTRAP, sizeof(double), sort_func);

	  /* IC basé sur les 2.5e et 97.5e percentiles */
	  printf("%s,%.8G,%.8G,%.8G,%.8G,%.8G,%.8G,%.8G,%+.8G,%+.8G,",
		 scenarios_list[scenarios_new_ranks[v]], c_cmp, inter_c_cmp,
		 delta_cmp, c_var, inter_c_var, delta_var, ICR,
		 ICR_bootstrap[(int)floor((double)BOOTSTRAP*0.025)]-ICR,
		 ICR_bootstrap[(int)ceil((double)BOOTSTRAP*0.975)]-ICR);

	  if (delta_var < 0)
	    {
	      printf("Dominé");
	      ++how_much_backward;
	    }
	  else
	    {
	      how_much_backward = 0;
	    }
	  printf("\n");
	}
    }
  return 0;
}

/*
 * Parcourt un fichier de configuration pour en extraire les spécifications
 * de l'usager.
 */
void  parse_configuration  (FILE * pConf, conf_args * to_fill,
			    int *vars_types, char **vars_list,
			    int vars_count)
{
  /* Cette fonction est le résultat de l'ajout successif de plusieurs
   * fonctionnalités. Elle pourrait être réécrite en utilisant GNU Bison.
   * Le comportement de la fonction, à l'opposé du code lui-même, devrait
   * être très aisé à comprendre, et il est de plus expliqué en détails
   * dans le manuel de ce programme.
   */

  /* Sauvegarder le point d'entrée */
  fpos_t orig;
  fgetpos (pConf, &orig);

  char  choices[][32] = { "proportions", "calculs (global)", "ICER",
			  "calculs (local)", "expressions booleennes",
			  "ne pas afficher", "calculs (conditionnel)",
			  {NUL} } ;
  char  line  [BUFFER_SIZE]; /* buffer */
  char  *pch   ; /* Pointeur du buffer */
  char  *pch_h ; /* Pointeur "helpeur" */
  char  *current = choices[7]; /* catégorie en cours de traitement */
  int   index  ;
  int   valid  ;

  int discrete_vars_count  = 0  ;

  int calcs_count          = 0  ;
  int calcs_vars_count          ;

  int loc_count            = 0  ;
  int total_loc_count      = 0  ;
  int higher_nb_vars       = 0  ;
  int loc_higher_nb        = 0  ;

  int ICR_vars_count       = 0  ;

  int bool_count           = 0  ;
  int bool_higher_nb       = 0  ;

  /* Compter le nombre d'éléments pour allocation des tableaux.
   * Un peu de traitement d'erreurs.
   */
  while (fgets(line, BUFFER_SIZE, pConf) != NULL)
    {
      pch = line;

      while ( isspace (*pch) && *pch != '\n')
	pch++;

      if (*pch == '\n')
	continue;

      if (*pch == '[')
	{
	  current = pch+1;
	  valid = 0;

	  for (index = 0; index < 7 /* magic number */; index++)
	    {
	      if (! strncmp (current, choices[index],
			     strlen(choices[index])))
		{
		  valid = 1;
		  break;
		}
	    }

	  if (! valid)
	    {
	      printf("Catégorie non reconnue: %s", line);
	      exit(1);
	    }
	}
      /* Si la ligne n'est pas un commentaire */
      else if (*pch != '#' && *pch != ';')
	{
	  if ( *current == NUL )
	    {
	      printf("Aucune catégorie n'a été initialisée: %s", line);
	      exit(1);
	    }
	  else
	    {
	      switch (index)
		{
		  /* Variables discrètes */
		case 0:
		  pch_h = pch;
		  do
		    pch_h++;
		  while ( ! isspace(*pch_h) );
		  *pch_h = NUL;

		  /* Changer le type */
		  vars_types [find_var_rank(pch, vars_count, vars_list)]
		    = DISCRETE;
		  discrete_vars_count++;
		  break;

		  /* Calculs entre variables (global) */
		case 1:
		  calcs_count++;
		  calcs_vars_count = 0;

		  pch = strchr (pch, '"');
		  while (pch != NULL)
		    {
		      pch++;
		      calcs_vars_count++;
		      pch = strchr (pch, '"');
		    }
		  calcs_vars_count >>= 1;

		  /* Trouver l'expression avec le plus grand nombre de
		   * variables: c'est l'allocation de mémoire qui sera
		   * utilisée pour tous les autres */
		  if (calcs_vars_count > higher_nb_vars)
		    higher_nb_vars = calcs_vars_count;

		  break;

		  /* ICR */
		case 2:
		  if (! strncmp (pch, "variable", 8))
		    {
		      ICR_vars_count++;
		    }
		  break;

		  /* Calculs entre variables (local et conditionnel)
		   * et expressions booléennes */
		case 3: case 4: case 6:
		  if (index == 3)
		    {
		      loc_count++;
		      total_loc_count++;
		    }

		  else if (index == 4)
		    bool_count++;

		  else if (! strncmp (pch, "condition", 9) )
		    continue;

		  else
		    total_loc_count++;

		  calcs_vars_count = 0;

		  pch = strchr (pch, '"');
		  while (pch != NULL)
		    {
		      pch++;
		      calcs_vars_count++;
		      pch = strchr (pch, '"');
		    }
		  calcs_vars_count >>= 1;

		  /* Encore une fois trouver le plus grand nombre de
		   * variables pour déterminer l'allocation de mém. */
		  if (index == 3 && calcs_vars_count > loc_higher_nb)
		    loc_higher_nb = calcs_vars_count;
		  else if (calcs_vars_count > bool_higher_nb)
		    bool_higher_nb = calcs_vars_count;

		  break;

		  /* Enlever l'affichage */
		case 5:
		  /* rien à faire */
		  break;
		}
	    }
	}
    }

  /* Ce qui est en malloc doit être transmis au struct contenant la
   * configuration, ce qui empêche d'en faire des variables locales. */

  char  **calcs_list    = (char**) malloc (sizeof(char*) * calcs_count);
  char  **calcs_labels  = (char**) malloc (sizeof(char*) * calcs_count);

  char  **loc_list      = (char**) malloc (sizeof(char*) * total_loc_count);
  char  **loc_labels    = (char**) malloc (sizeof(char*) * total_loc_count);

  char  **bool_labels   = (char**) malloc (sizeof(char*) * bool_count);

  int   **calcs_relative_ranks = (int**) malloc (sizeof(int*)
						 * calcs_count);
  int   **calcs_vars_types     = (int**) malloc (sizeof(int*)
						 * calcs_count);
  for (int i = 0; i < calcs_count; ++i)
    {
      calcs_relative_ranks[i] = (int*) malloc (sizeof(int)
					       * higher_nb_vars);
      calcs_vars_types[i]     = (int*) malloc (sizeof(int)
					       * higher_nb_vars);
      calcs_labels[i]         = NULL;
    }

  int   *loc_vars_count  = (int*)  malloc (sizeof(int)  * total_loc_count);
  int   **loc_vars_ranks = (int**) malloc (sizeof(int*) * total_loc_count);
  for (int i = 0; i < total_loc_count; ++i)
    {
      loc_vars_ranks[i] = (int*) malloc (sizeof(int) * loc_higher_nb);
      loc_labels[i]     = NULL;
    }

  int  *cond_vars_rank  = (int*) malloc (sizeof(int) * total_loc_count
					 - loc_count);

  int   *bool_vars_count  = (int*)  malloc (sizeof(int)   * bool_count);
  int   **bool_vars_ranks = (int**) malloc (sizeof(int*)  * bool_count);
  int   **bool_op_list    = (int**) malloc (sizeof(int*)  * bool_count);
  int   **comp_op_list    = (int**) malloc (sizeof(int*)  * bool_count);
  char  ***data_comp_list = (char***) malloc (sizeof(char**) * bool_count);

  for (int i = 0; i < bool_count; ++i)
    {
      bool_vars_ranks[i] = (int*) malloc (sizeof(int) * bool_higher_nb);
      bool_op_list[i]    = (int*) malloc (sizeof(int) * bool_higher_nb);
      comp_op_list[i]    = (int*) malloc (sizeof(int) * bool_higher_nb);
      data_comp_list[i]  = (char**) malloc(sizeof(char*) * bool_higher_nb);
    }

  int   *ICR_vars_ranks   = (int*) malloc (sizeof(int) * ICR_vars_count);
  int   *ICR_vars_types   = (int*) malloc (sizeof(int) * ICR_vars_count);
  int   *ICR_vars_inv     = (int*) malloc (sizeof(int) * ICR_vars_count);

  int   *no_show          = (int*) calloc (vars_count + bool_count
					   + total_loc_count + calcs_count,
					   sizeof(int) );

  int   ICR_cmp_rank  = -1;
  int   ICR_cmp_type;
  int   inv_ready     =  0;
  int   ICR_rank      =  0;

  int   calcs_rank    =  0;
  int   loc_rank  ;
  int   bool_rank     =  0;

  int   bool_ready    =  0;
  int   cond_rank     =  loc_count;
  int   loc_norm_rank =  0;

  int   var_abs_rank;
  int   var_relative_rank;
  int   var_type;

  int bool_op_count;
  int comp_op_count;

  char to_save;
  char *equ_start;

  fsetpos (pConf, &orig); /* Revenir au début du fichier */

  /* Véritable traitement cette fois: remplir les tableaux alloués */
  while (fgets(line, BUFFER_SIZE, pConf) != NULL)
    {
      pch = line;

      while ( isspace (*pch) && *pch != '\n')
	pch++;

      if (*pch == '\n')
	continue;

      /* Moins de traitement d'erreurs car on s'est déjà assuré de la
       * validité de toutes les catégories */
      if (*pch == '[')
	{
	  current = pch+1;
	  for (index = 0; index < 7 /* magic number */; index++)
	    {
	      if (! strncmp (current, choices[index],
			     strlen(choices[index])))
		break;
	    }
	}

      else if (*pch != '#' && *pch != ';')
	{
	  switch (index)
	    {
	      /* variables discrètes */
	    case 0:
	      /* rien à faire */
	      break;

	      /* Calculs entre variables (global) */
	    case 1:
	      /* vérifie si une étiquette existe */
	      pch_h = strchr (pch, '=');
	      if (pch_h != NULL)
		{
		  pch_h = pch;
		  do
		    pch_h++;
		  while ( !isspace (*pch_h) && *pch_h != '=');

		  to_save = *pch_h;
		  *pch_h = NUL;

		  calcs_labels[calcs_rank] = (char*) malloc (strlen(pch)
							     + 1);
		  strcpy (calcs_labels[calcs_rank], pch);

		  *pch_h = to_save;

		  pch = pch_h;

		  do
		    pch++;
		  while ( ( isspace (*pch) || *pch == '=' )
			  && *pch != '\n' );
		}

	      /* Lecture de l'équation */
	      equ_start = pch;

	      calcs_vars_count = 0;
	      while (*pch != '\n')
		{
		  if (*pch == '"')
		    {
		      pch_h = strchr(pch+1, '"');
		      if (pch_h == NULL)
			{
			  printf("Une variable n'a pas été refermée par un \
guillemet dans les calculs: %s", line) ;
			  exit(1);
			}
		      *pch_h = NUL;

		      find_var_type_and_rank(++pch, vars_count, vars_list,
					     vars_types, calcs_count,
					     calcs_labels, total_loc_count,
					     loc_labels, bool_count,
					     bool_labels, &var_type,
					     &var_abs_rank);

		      if (var_type == DISCRETE)
			{
			  printf("Impossible d'utiliser la variable '%s'. \
Type invalide.\n", pch);
			  exit(1);
			}

		      /* Si variable de type standard, il faut trouver
		       * sont rang relatif, c.a.d. par rapport à son propre
		       * type de variable. */
		      if (var_type < DISCRETE)
			{
			  var_relative_rank = 0;
			  for (int i = 0; i < var_abs_rank; i++)
			    if (vars_types[i] == vars_types[var_abs_rank])
			      var_relative_rank++;

			  var_abs_rank = var_relative_rank;
			}

		      calcs_relative_ranks[calcs_rank][calcs_vars_count]
			= var_abs_rank;
		      calcs_vars_types[calcs_rank][calcs_vars_count++]
			= var_type;

		      *pch_h = '"';
		      pch = pch_h + 1;
		    }
		  else if (! isspace (*pch) && ! isdigit (*pch)
			   && *pch != '.')
		    {
		      check_delim_exist(*pch);
		      ++pch;
		    }
		  else
		    ++pch;
		}

	      /* On copie le calcul lui-même après avoir vérifié la
	       * validité de toutes les variables et enregistrer leur
	       * rang et type. */
	      calcs_list[calcs_rank] = (char*) malloc(strlen(equ_start));
	      strncpy(calcs_list[calcs_rank], equ_start, strlen(equ_start));
	      calcs_list[calcs_rank++][strlen(equ_start)-1] = NUL;
	      break;

	      /* ICR: Trouver le rang et type des variables et du
	       * comparateur. Traiter aussi l'inversion. */
	    case 2:
	      /* Variable */
              pch = strchr (line, '=');
              if (pch == NULL)
                {
                  printf("Aucun symbole '=': %s", line);
                  exit(1);
                }

	      else if (! strncmp (line, "variable", 8))
		{
		  if (inv_ready)
		    {
		      printf("Il faut préciser l'inversion de la variable \
'%s' avant de passer à une autre variable.\n",
			     vars_list[ICR_vars_ranks[ICR_rank]]);
                      exit(1);
		    }

		  do
		    pch++;
		  while ( isspace (*pch) );

		  pch_h = pch;
		  do
		    pch_h++;
		  while ( ! isspace(*pch_h) );
		  *pch_h = NUL;

		  find_var_type_and_rank(pch, vars_count, vars_list,
					 vars_types, calcs_count,
					 calcs_labels, total_loc_count,
					 loc_labels, bool_count,
					 bool_labels,
					 ICR_vars_types+ICR_rank,
					 ICR_vars_ranks+ICR_rank);

		  inv_ready = 1;
		}
	      /* Inversion ou non d'une variable */
	      else if (! strncmp(line, "inv", 3))
		{
		  if (! inv_ready)
		    {
		      puts("Il faut préciser une variable avant de \
spécifier une inversion.");
		      exit(1);
                    }

		  do
		    pch++;
		  while ( isspace (*pch) );

		  pch_h = pch;
		  do
		    pch_h++;
		  while ( isdigit (*pch_h) );
		  *pch_h = NUL;

		  ICR_vars_inv[ICR_rank++] = atoi(pch);
		  inv_ready = 0;
		}
	      /* Comparateur */
	      else if (! strncmp(line, "comparateur", 11))
		{
		  if (ICR_cmp_rank >= 0)
		    {
		      printf("Un comparateur existe déjà: '%s'\n",
			     vars_list[ICR_cmp_rank]) ;
		      exit(1);
		    }

		  do
		    pch++;
		  while ( isspace (*pch) );

		  pch_h = pch;
		  do
		    pch_h++;
		  while ( ! isspace (*pch_h) );
		  *pch_h = NUL;

		  find_var_type_and_rank(pch, vars_count, vars_list,
					 vars_types, calcs_count,
					 calcs_labels, total_loc_count,
					 loc_labels, bool_count,
					 bool_labels, &ICR_cmp_type,
					 &ICR_cmp_rank);
		}
	      else
		{
		  printf("Argument non reconnu: %s", line);
		  exit(1);
		}
	      break;

	      /* Calculs entre variables (local et conditionnel):
	       * similaire aux calculs globaux, mais plus d'information
	       * est nécessaire. */
	    case 3: case 6:
	      if (index == 6 && !strncmp (pch, "condition", 9))
		{
		  if (bool_ready)
		    {
		      pch = strchr (line, '=');
                      if (pch == NULL)
                        {
                          printf("Aucun symbole '=': %s", line);
                          exit(1);
                        }

		      do
			pch++;
		      while ( isspace (*pch) );

		      pch_h = pch;
		      do
			pch_h++;
		      while ( !isspace (*pch_h) );
		      *pch_h = NUL;

		      find_var_type_and_rank(pch, vars_count, vars_list,
					     vars_types, calcs_count,
					     calcs_labels, total_loc_count,
					     loc_labels, bool_count,
					     bool_labels, &var_type,
					     &var_relative_rank);

		      if (! (var_type == CUSTOM_BOOLEAN
			     || var_type == BOOLEAN) )
			{
			  printf("Type invalide: %s\n", pch);
			  exit(1);
			}

		      cond_vars_rank[cond_rank - loc_count] =
			var_relative_rank;
		      if (var_type == CUSTOM_BOOLEAN)
			cond_vars_rank[cond_rank - loc_count] += vars_count
			  + loc_count;

		      bool_ready = 0;
		      cond_rank++;
		    }
		  else
		    {
		      printf("Il faut préciser un calcul avant de spécifier\
 une condition: %s", line);
		      exit(1);
		    }
		}
	      else if (bool_ready)
		{
		  printf("Il faut préciser la condition du calcul '%s' \
avant de passer à un calcul suivant.\n", loc_list[cond_rank]);
		  exit(1);
		}
	      else
		{
		  loc_rank = index == 3 ? loc_norm_rank : cond_rank ;
		  /* vérifie si une étiquette existe */
		  pch_h = strchr (pch, '=');
		  if (pch_h != NULL)
		    {
		      pch_h = pch;
		      do
			pch_h++;
		      while ( !isspace (*pch_h) && *pch_h != '=');

		      to_save = *pch_h;
		      *pch_h = NUL;

		      loc_labels[loc_rank] = (char*) malloc (strlen(pch)
							     + 1);
		      strcpy (loc_labels[loc_rank], pch);

		      *pch_h = to_save;

		      pch = pch_h;

		      do
			pch++;
		      while ( ( isspace (*pch) || *pch == '=' )
			      && *pch != '\n' ) ;
		    }

		  /* Traitement de l'équation ici */
		  equ_start = pch;
		  calcs_vars_count = 0;
		  while (*pch != '\n')
		    {
		      if (*pch == '"')
			{
			  pch_h = strchr(pch+1, '"');
			  if (pch_h == NULL)
			    {
			      printf("Une variable n'a pas été refermée par\
 un guillemet dans les calculs: %s", line) ;
			      exit(1);
			    }
			  *pch_h = NUL;
			  find_var_type_and_rank(++pch, vars_count,
						 vars_list, vars_types,
						 calcs_count, calcs_labels,
						 total_loc_count,
						 loc_labels, bool_count,
						 bool_labels, &var_type,
						 &var_relative_rank);

			  /* Les types booléens ne sont pas acceptés dans
			   * les calculs. */
			  if (var_type == BOOLEAN
			      || var_type == CUSTOM_BOOLEAN
			      || var_type == GLOB_CALC)
			    {
			      printf("Impossible d'utiliser la variable \
'%s'. Type invalide.\n", pch) ;
			      exit(1);
			    }

			  loc_vars_ranks[loc_rank][calcs_vars_count++]
			    = var_type == LOC_CALC ? var_relative_rank +
			    vars_count : var_relative_rank ;

			  *pch_h = '"';
			  pch = pch_h + 1;
			}
		      else if (! isspace (*pch) && ! isdigit (*pch))
			{
			  check_delim_exist(*pch);
			  ++pch;
			}
		      else
			++pch;
		    }

		  if (! calcs_vars_count)
		    {
		      printf("Calcul invalide! Il faut utiliser au moins \
une variable: %s", line) ;
		      exit(1);
		    }

		  /* Tout est valide: copie de l'équation en mémoire */
		  loc_vars_count[loc_rank] = calcs_vars_count;
		  loc_list[loc_rank] = (char*) malloc(strlen(equ_start));
		  strncpy(loc_list[loc_rank], equ_start, strlen(equ_start));
		  loc_list[loc_rank][strlen(equ_start)-1] = NUL;

		  if (index == 3)
		    loc_norm_rank++;
		  else
		    bool_ready = 1;
		}
	      break;

	      /* Expressions booléennes */
	    case 4:
	      /* vérifie si une étiquette existe */
	      pch_h = strchr (pch, '=');
	      if (pch_h != NULL && *(pch_h+1) != '=' && *(pch_h+1) != '!')
		{
		  pch_h = pch;
		  do
		    pch_h++;
		  while ( !isspace (*pch_h) && *pch_h != '=');

		  to_save = *pch_h;
		  *pch_h = NUL;

		  bool_labels[bool_rank] = (char*) malloc (strlen(pch) + 1);
		  strcpy (bool_labels[bool_rank], pch);

		  *pch_h = to_save;

		  pch = pch_h;

		  do
		    pch++;
		  while ( ( isspace (*pch) || *pch == '=' )
			  && *pch != '\n' );
		}
	      else
		{
		  printf("Il faut donner une étiquette à l'expression \
booléenne: %s", line) ;
		  exit(1);
		}

	      calcs_vars_count = 0;
	      bool_op_count    = 0;
	      comp_op_count    = 0;

	      while (*pch != '\n')
		{
		  if (*pch == '"')
		    {
		      pch_h = strchr(pch+1, '"');
		      if (pch_h == NULL)
			{
			  printf("Une variable n'a pas été refermée par un \
guillemet dans les expressions booléennes: %s",
				 line) ;
			  exit(1);
			}
		      *pch_h = NUL;

		      find_var_type_and_rank(++pch, vars_count, vars_list,
					     vars_types, calcs_count,
					     calcs_labels, total_loc_count,
					     loc_labels, bool_count,
					     bool_labels, &var_type,
					     &var_relative_rank);

		      /* Les calculs globaux ne sont calculés qu'après le
		       * parsing... Utiliser les calculs locaux. */
		      if (var_type == GLOB_CALC)
			{
			  printf("Impossible d'utiliser la variable '%s'. \
Type invalide.\n",
				 pch);
			  exit(1);
			}

		      if (var_type <= DISCRETE)
			bool_vars_ranks[bool_rank][calcs_vars_count]
			  = var_relative_rank;
		      else
			bool_vars_ranks[bool_rank][calcs_vars_count]
			  = var_type == LOC_CALC ? vars_count
			  + var_relative_rank : vars_count + loc_count
			  + var_relative_rank;

		      *pch_h = '"';
		      pch = pch_h + 1;
		    }
		  else if (! isspace (*pch) && ! isdigit (*pch))
		    {
		      /* Traite les opérateurs trouvés */
		      if (*pch == '&' || *pch == '|')
			{
			  if ( (*pch == '&' && *(pch+1) != '&')
			      || (*pch == '|' && *(pch+1) != '|')
			      || calcs_vars_count-1 != bool_op_count )
			    {
			      printf("Expression invalide: %s", line);
			      exit(1);
			    }
			  bool_op_list[bool_rank][calcs_vars_count - 1]
			    = *pch == '&' ? AND : OR ;
			  bool_op_count++;
			  pch += 2;
			}
		      else if (*pch == '<' || *pch == '>' || *pch == '='
			       || *pch == '!')
			{
			  if (calcs_vars_count != comp_op_count)
			    {
			      printf("Expression invalide: %s", line);
			      exit(1);
			    }
			  if (*pch == '=' || *pch == '!')
			    {
			      if ( *(pch+1) != '=')
				{
				  printf("Expression invalide: %s", line);
				  exit(1);
				}

			      comp_op_list[bool_rank][calcs_vars_count]
				= *pch == '=' ? EQ : NE ;
			      ++pch;
			    }
			  else if (*(pch+1) == '=')
			    {
			      comp_op_list[bool_rank][calcs_vars_count]
				= *pch == '<' ? LE : GE ;
			      ++pch;
			    }
			  else
			    {
			      comp_op_list[bool_rank][calcs_vars_count]
				= *pch == '<' ? LT : GT ;
			    }

			  comp_op_count++;
			  do
			    ++pch;
			  while (*pch != '\n' && isspace (*pch) );

			  if (*pch == '\n')
			    {
			      printf("Expression invalide: %s", line);
			      exit(1);
			    }

			  pch_h = pch;
			  do
			    ++pch_h;
			  while (! isspace (*pch_h) && *pch_h != '&'
				 && *pch_h != '|') ;

			  to_save = *pch_h;
			  *pch_h = NUL;

			  data_comp_list[bool_rank][calcs_vars_count]
			    = (char*) malloc ( strlen(pch) + 1);

			  strcpy (data_comp_list[bool_rank]
				  [calcs_vars_count++], pch) ;

			  *pch_h = to_save;
			  pch = pch_h;
			}
		      else
			{
			  printf("Expression invalide: %s", line);
			  exit(1);
			}
		    }
		  else
		    ++pch;
		}

	      if (! calcs_vars_count)
		{
		  printf("Expression invalide! Il faut utiliser au moins \
une variable: %s", line) ;
		  exit(1);
		}

	      bool_vars_count[bool_rank++] = calcs_vars_count;
	      break;

	      /* Enlever l'affichage */
	    case 5:
	      pch_h = pch;
	      do
		++pch_h;
	      while ( ! isspace (*pch_h) );

	      *pch_h = NUL;
	      find_var_type_and_rank(pch, vars_count, vars_list, vars_types,
				     calcs_count, calcs_labels,
				     total_loc_count, loc_labels,
				     bool_count, bool_labels, &var_type,
				     &var_relative_rank) ;

	      switch (var_type)
		{
		case BOOLEAN: case ACCUMUL: case DISCRETE:
		  no_show [var_relative_rank] = 1;
		  break;

		case CUSTOM_BOOLEAN:
		  no_show [vars_count + var_relative_rank] = 1;
		  break;

		case LOC_CALC:
		  no_show [vars_count + bool_count + var_relative_rank] = 1;
		  break;

		case GLOB_CALC:
		  no_show [vars_count + bool_count + total_loc_count
			   + var_relative_rank] = 1;
		  break;
		}
	      break;
	    }
	}
    }

  /* Remplir le struct */
  to_fill->calcs_list     = calcs_list     ;
  to_fill->calcs_labels   = calcs_labels   ;
  to_fill->loc_list       = loc_list       ;
  to_fill->loc_labels     = loc_labels     ;

  to_fill->calcs_relative_ranks = calcs_relative_ranks ;
  to_fill->calcs_vars_types     = calcs_vars_types     ;
  to_fill->loc_vars_ranks       = loc_vars_ranks       ;
  to_fill->loc_vars_count       = loc_vars_count       ;

  to_fill->bool_vars_count      = bool_vars_count  ;
  to_fill->bool_vars_ranks      = bool_vars_ranks  ;
  to_fill->bool_op_list         = bool_op_list     ;
  to_fill->comp_op_list         = comp_op_list     ;
  to_fill->data_comp_list       = data_comp_list   ;
  to_fill->bool_labels          = bool_labels      ;

  to_fill->cond_vars_rank       = cond_vars_rank   ;
  to_fill->no_show              = no_show          ;

  to_fill->ICR_vars_ranks = ICR_vars_ranks ;
  to_fill->ICR_cmp_rank   = ICR_cmp_rank   ;
  to_fill->ICR_cmp_type   = ICR_cmp_type   ;
  to_fill->ICR_vars_types = ICR_vars_types ;
  to_fill->ICR_vars_inv = ICR_vars_inv ;

  to_fill->discrete_vars_count = discrete_vars_count ;
  to_fill->calcs_count         = calcs_count         ;
  to_fill->loc_count           = loc_count           ;
  to_fill->total_loc_count     = total_loc_count     ;
  to_fill->bool_count          = bool_count          ;
  to_fill->ICR_vars_count      = ICR_vars_count      ;

  if (ICR_vars_count && ICR_cmp_rank < 0)
    {
      puts("Des variables d\'ICR existent, mais aucun comparateur n'a été \
donné.");
      exit(1);
    }

  else if (inv_ready)
    {
      printf("Il faut préciser l'inversion de la variable '%s'\n",
	     vars_list[ICR_vars_ranks[ICR_rank]]) ;
      exit (1);
    }

  else if (bool_ready)
    {
      printf("Il faut préciser la condition du calcul: %s\n",
	     loc_list[cond_rank]);
      exit(1);
    }

  else
    return;
}

/*
 * Retourne le rang de la variable 'var'. (variables standards seulement)
 */
int find_var_rank (char *var, int vars_count, char **vars_list)
{
  for (int i = 0; i < vars_count; ++i)
    {
      if (! strcmp (vars_list[i], var))
	{
	  return i;
	}
    }
  printf("'%s' n'est pas une variable reconnue. ", var);
  printf("Voici une liste des variables disponibles:\n");
  for (int i = 0; i < vars_count; i++)
    printf(" - %s\n", vars_list[i]);
  exit(1);
}

/*
 * Trouve le type et le rang de la variable 'var'. (toutes les variables
 * traitées jusqu'à l'appel de la fonction)
 */
void find_var_type_and_rank (char *var, int vars_count, char **vars_list,
			     int *vars_types, int calcs_count,
			     char **calcs_labels, int loc_count,
			     char **loc_labels, int bool_count,
			     char **bool_labels, int *type, int *rank)
{
  for (int i = 0; i < vars_count; ++i)
    {
      if (! strcmp (vars_list[i], var))
	{
	  *type = vars_types[i];
	  *rank = i;
	  return;
	}
    }
  for (int i = 0; i < calcs_count; ++i)
    {
      if (calcs_labels[i] != NULL)
	if (! strcmp (calcs_labels[i], var))
	  {
	    *type = GLOB_CALC;
	    *rank = i;
	    return;
	  }
    }
  for (int i = 0; i < loc_count; ++i)
    {
      if (loc_labels[i] != NULL)
	if (! strcmp (loc_labels[i], var))
	  {
	    *type = LOC_CALC;
	    *rank = i;
	    return;
	  }
    }
  for (int i = 0; i < bool_count; ++i)
    {
      if (! strcmp (bool_labels[i], var))
	{
	  *type = CUSTOM_BOOLEAN;
	  *rank = i;
	  return;
	}
    }
  printf("'%s'n'est pas une variable reconnue. ", var);
  printf("Voici une liste des variables disponibles:\n");
  int i;
  for (i = 0; i < vars_count; i++)
    printf(" - %s\n", vars_list[i]);
  for (i = 0; i < calcs_count; i++)
    if (calcs_labels[i] != NULL)
      printf(" - %s\n", calcs_labels[i]);
  for (i = 0; i < loc_count; i++)
    if (loc_labels[i] != NULL)
      printf(" - %s\n", loc_labels[i]);
  for (i = 0; i < bool_count; i++)
    printf(" - %s\n", bool_labels[i]);
  exit(1);
}

/*
 * Vérifie si le délimiteur existe. (pour les calculs)
 */
void check_delim_exist (char delim)
{
  for (int i = 0; i < 7; i++)
    {
      /* delims est un tableau à portée globale défini dans eval.h: le fait
       * qu'il soit global provient du code original */
      if (delim == delims[i])
	return;
    }
  printf("Le délimiteur '%c' dans les calculs n'est pas valide.\n", delim);
  exit(1);
}

/*
 * Parcourt les fichiers CSV.
 */
void * parse_csv (void *ptr)
{
  /* Transformer un pointeur void en pointeur de struct */
  thread_args *struct_Ptr = (thread_args *) ptr;

  int str_length  = strlen (struct_Ptr->path) + strlen (struct_Ptr->name);
  char file_path [BUFFER_SIZE] ;
  strcpy (file_path, struct_Ptr->path);
  strcat (file_path, struct_Ptr->name);

  char line [BUFFER_SIZE] ;
  char iter_buffer [32] ;
  char *dump ; /* Pour rendre une fonction thread-safe */

  int i, v, k, p ;

  int acc_vars_rank  ;
  int bool_vars_rank ;
  int dis_vars_rank  ;

  int current_c_bool_state;
  int previous_c_bool_state;

  /* offset pour variables accumulatrices */
  int offset_acc = ( ( struct_Ptr->iters_count * struct_Ptr->acc_vars_count
		       * struct_Ptr->num_scen )
		     + ( struct_Ptr->lower_lim
			 * struct_Ptr->acc_vars_count ) );

  /* offset pour variables booléennes */
  int offset_bo = ( ( struct_Ptr->iters_count * struct_Ptr->bool_vars_count
		    * struct_Ptr->num_scen )
		   + ( struct_Ptr->lower_lim
		       * struct_Ptr->bool_vars_count ) );

  /* offset pour variables discrètes */
  int offset_dis = ( ( struct_Ptr->iters_count
		       * struct_Ptr->discrete_vars_count
		       * struct_Ptr->num_scen )
		     + ( struct_Ptr->lower_lim
			 * struct_Ptr->discrete_vars_count ) );

  /* offset pour calculs */
  int offset_loc = ((struct_Ptr->iters_count * struct_Ptr->total_loc_count
		       * struct_Ptr->num_scen )
		     + ( struct_Ptr->lower_lim
			 * struct_Ptr->total_loc_count) );

  /* offset pour booléens personnalisés */
  int offset_c_bo = ( ( struct_Ptr->iters_count * struct_Ptr->c_bool_count
			* struct_Ptr->num_scen )
		      + ( struct_Ptr->lower_lim
			  * struct_Ptr->c_bool_count ) );

  /* On initialise sa part des tableaux de résultats à 0 */
  memset ( struct_Ptr->acc_results + offset_acc, 0,
	   sizeof(double) * ( ( struct_Ptr->upper_lim
				- struct_Ptr->lower_lim )
			      * struct_Ptr->acc_vars_count ) );

  memset ( struct_Ptr->bool_results + offset_bo, 0,
	   sizeof(unsigned int) * ( ( struct_Ptr->upper_lim
				      - struct_Ptr->lower_lim )
				    * struct_Ptr->bool_vars_count ) );

  memset ( struct_Ptr->loc_results + offset_loc, 0,
	   sizeof(double) * ( (struct_Ptr->upper_lim
			       - struct_Ptr->lower_lim)
			      * struct_Ptr->total_loc_count ) );

  memset ( struct_Ptr->c_bool_results + offset_c_bo, 0,
	   sizeof(unsigned int) * ( (struct_Ptr->upper_lim
				     - struct_Ptr->lower_lim)
				    * struct_Ptr->c_bool_count ) );

  /* Évite une multiplication dans la boucle for */
  offset_acc  -= struct_Ptr->acc_vars_count      ;
  offset_bo   -= struct_Ptr->bool_vars_count     ;
  offset_dis  -= struct_Ptr->discrete_vars_count ;
  offset_loc  -= struct_Ptr->total_loc_count     ;
  offset_c_bo -= struct_Ptr->c_bool_count        ;

  string value;
  double num_value;

  /* Une cache qui contient les informations sur le dernier individu. */
  last_value *cache = (last_value*) malloc
    (sizeof(last_value)
     * (struct_Ptr->vars_count + struct_Ptr->total_loc_count
	+ struct_Ptr->c_bool_count));

  char tokenizer [BUFFER_SIZE]; /* pour recevoir une copie qui sera
				 * mise en morceaux par strtok_r */
  char buffer    [BUFFER_SIZE];

  char *pch, *dump_calc;

  int return_check;
  eval evaluator; /* Permet les calculs. */

  /* Évite un avertissement de conversion dépréciée */
  char true_s[] = "true";
  char false_s[] = "false";

  char   *parsing;
  gzFile p_file;

  for (i = struct_Ptr->lower_lim; i < struct_Ptr->upper_lim; ++i)
    {
      /* Ouvrir le fichier de la simulation à analyser */
      sprintf (iter_buffer, "/%d_Output.gz", i);
      strcpy (file_path + str_length, iter_buffer);

      p_file = gzopen (file_path, "rb");

      if (p_file == Z_NULL)
	{
	  printf("Mauvais chemin ou fichier inexistant: %s\n", file_path);
	  exit(1);
	}

      gzgets (p_file, line, BUFFER_SIZE);

      /* Mettre à jour l'offset */
      offset_acc  += struct_Ptr->acc_vars_count      ;
      offset_bo   += struct_Ptr->bool_vars_count     ;
      offset_dis  += struct_Ptr->discrete_vars_count ;
      offset_loc  += struct_Ptr->total_loc_count     ;
      offset_c_bo += struct_Ptr->c_bool_count        ;

      /* Parsing selon la population et la colonne (fichier CSV) */
      for (v = 0; v < struct_Ptr->pop; ++v)
	{
	  gzgets(p_file, line, BUFFER_SIZE);
	  strtok_r(line, ",", &dump);

	  acc_vars_rank  = 0;
	  bool_vars_rank = 0;
	  dis_vars_rank  = 0;

	  /* Les variables standards */
	  for (k = 0; k < struct_Ptr->vars_count; ++k)
	    {
	      parsing = strtok_r(NULL, ",", &dump);

	      if (parsing == NULL)
		{
		  printf("Inexistant! Fichier: %s, Ligne: %d, Var #: %d\n",
			 file_path, v, k);
		  exit(1);
		}

	      /* Enregistrer les résultats dans le tableau approprié
	       * et la cache */
	      switch ( struct_Ptr->vars_types [k] )
		{
		case BOOLEAN:
		  if (! strncmp(parsing, "true", 4))
		      ++struct_Ptr->bool_results[offset_bo+bool_vars_rank];

		  cache[k].string_value = parsing ;
		  ++bool_vars_rank;
		  break;

		case ACCUMUL:
		  num_value = atof (parsing) ;
		  struct_Ptr->acc_results [offset_acc + acc_vars_rank++]
		    += num_value ;
		  cache[k].num_value = num_value ;
		  break;

		case DISCRETE:
		  value = string (parsing);
		  if (value[value.size() - 1] == '\n')
		    value[value.size() - 1] = NUL;
		  ++struct_Ptr->discrete_results[offset_dis
						 + dis_vars_rank++][value];
		  cache[k].string_value = parsing ;
		  break;
		}
	    }

	  /* Calculs locaux */
	  for (k = 0; k < struct_Ptr->loc_count; ++k)
	    {
	      strcpy (tokenizer, struct_Ptr->loc_list[k]);

	      if (tokenizer[0] != "\""[0])
		{
		  pch = strtok_r (tokenizer, "\"", &dump_calc);
		  strcpy (buffer, pch);
		  pch = strtok_r (NULL, "\"", &dump_calc);
		}
	      else
		{
		  strcpy (buffer, "");
		  pch = strtok_r(tokenizer, "\"", &dump_calc);
		  pch++;
		}

	      if (struct_Ptr->vars_types[struct_Ptr->loc_vars_ranks[k][0]]
		  == DISCRETE)
		cache[struct_Ptr->loc_vars_ranks[k][0]].num_value
		  = atof(cache
			 [struct_Ptr->loc_vars_ranks[k][0]].string_value);

	      sprintf (buffer, "%s%f", buffer,
		       cache[struct_Ptr->loc_vars_ranks[k][0]].num_value);

	      for (p = 1; p < struct_Ptr->loc_vars_count[k]; ++p)
		{
		  pch = strtok_r (NULL, "\"", &dump_calc);
		  strcat (buffer, pch);
		  pch = strtok_r (NULL, "\"", &dump_calc);

		  sprintf (buffer, "%s%f", buffer,
			   cache[struct_Ptr->loc_vars_ranks[k][p]]
			   .num_value);
		}

	      pch = strtok_r (NULL, "\"", &dump_calc);

	      if (pch != NULL)
		strcat(buffer, pch);

	      /* appel d'un équivalent C du 'eval' pythonesqe */
	      return_check = evaluator.evaluate (buffer, &num_value);

	      if (return_check == -2)
		{
		  printf("Erreur de champ (range) dans les calculs locaux: \
%s, calcul: %s, scénario: %s, iteration %d\n",
			 struct_Ptr->loc_list[k], buffer,
			 struct_Ptr->name, i);
		  exit(1);
		}

	      cache [struct_Ptr->vars_count + k].num_value = num_value;
	      struct_Ptr->loc_results [offset_loc + k] += num_value ;
	    }

	  /* Expressions booléennes */
	  for (k = 0; k < struct_Ptr->c_bool_count; ++k)
	    {
	      for (p = 0; p < struct_Ptr->c_bool_vars_count[k]; ++p)
		{
		  previous_c_bool_state = current_c_bool_state;

		  /* Traiter l'opérateur rencontré */
		  switch (struct_Ptr->comp_op_list[k][p])
		    {
		    case EQ:
		      if ((struct_Ptr->c_bool_vars_ranks[k][p]
			   >= struct_Ptr->vars_count
			   && struct_Ptr->c_bool_vars_ranks[k][p]
			   < struct_Ptr->vars_count + struct_Ptr->loc_count)
			  || struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]] == ACCUMUL)
			{
			  if ( cache[struct_Ptr->c_bool_vars_ranks[k][p]]
			       .num_value == atof
			       (struct_Ptr->data_comp_list[k][p]) )
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (! strcmp (cache [struct_Ptr->
					       c_bool_vars_ranks[k][p]]
					.string_value, struct_Ptr->
					data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;

		    case NE:
		      if ((struct_Ptr->c_bool_vars_ranks[k][p]
			   >= struct_Ptr->vars_count
			   && struct_Ptr->c_bool_vars_ranks[k][p]
			   < struct_Ptr->vars_count + struct_Ptr->loc_count)
			  || struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]] == ACCUMUL)
			{
			  if (cache [struct_Ptr->c_bool_vars_ranks[k][p]]
			      .num_value != atof
			      (struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (strcmp (cache [struct_Ptr->c_bool_vars_ranks
					     [k][p]].string_value,
				      struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;

		    case GT:
		      if (struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]]
			  == DISCRETE)
			{
			  if (atof(cache[struct_Ptr->c_bool_vars_ranks
					 [k][p]].string_value)
			      > atof(struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (cache [struct_Ptr->c_bool_vars_ranks[k][p]]
			      .num_value > atof
			      (struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;

		    case GE:
		      if (struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]] == DISCRETE)
			{
			  if (atof(cache [struct_Ptr->c_bool_vars_ranks
					  [k][p]].string_value)
			      >= atof(struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (cache [struct_Ptr->c_bool_vars_ranks[k][p]]
			      .num_value >= atof(struct_Ptr->data_comp_list
						 [k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;

		    case LT:
		      if (struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]] == DISCRETE)
			{
			  if (atof(cache [struct_Ptr->c_bool_vars_ranks
					  [k][p]].string_value)
			      < atof(struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (cache [struct_Ptr->c_bool_vars_ranks[k][p]]
			      .num_value < atof
			      (struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;

		    case LE:
		      if (struct_Ptr->vars_types
			  [struct_Ptr->c_bool_vars_ranks[k][p]] == DISCRETE)
			{
			  if (atof(cache [struct_Ptr->c_bool_vars_ranks
					  [k][p]].string_value)
			      <= atof(struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      else
			{
			  if (cache [struct_Ptr->c_bool_vars_ranks[k][p]]
			      .num_value <= atof
			      (struct_Ptr->data_comp_list[k][p]))
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			}
		      break;
		    }
		  /* Vérifie si l'on est rendu à un délimiteur de type
		     && ou || */
		  if (p & 1)
		    {
		      /* Mettre à jour la valeur actuelle en prenant en
		       * en compte la valeur actuelle et précédente */
		      switch (struct_Ptr->bool_op_list[k][p >> 1])
			{
			case AND:
			  if (previous_c_bool_state && current_c_bool_state)
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			  break;

			case OR:
			  if (previous_c_bool_state || current_c_bool_state)
			    current_c_bool_state = 1;
			  else
			    current_c_bool_state = 0;
			  break;
			}
		    }
		}
	      if (current_c_bool_state)
		{
		  ++struct_Ptr->c_bool_results [offset_c_bo + k] ;
		  cache[struct_Ptr->vars_count+struct_Ptr->loc_count+k]
		    .string_value = true_s;
		}
	      else
		cache[struct_Ptr->vars_count+struct_Ptr->loc_count+k]
		  .string_value = false_s;
	    }

	  /* Calculs avec condition */
	  for (k = struct_Ptr->loc_count; k < struct_Ptr->total_loc_count;
	       ++k)
	    {
	      if ( strncmp(cache[struct_Ptr->cond_vars_rank
				[k - struct_Ptr->loc_count]]
			   .string_value, "true", 4) )
		{
		  /* La condition est fausse: rien à faire sauf mettre la
		   * cache à jour. */
		  cache [struct_Ptr->vars_count + struct_Ptr->c_bool_count
			 + k].num_value = 0;
		}

	      /* Simple copier coller des calculs locaux. (on ne peut pas
	       * les jumeler parce que ceux avec conditions doivent suivre
	       * les expressions booléennes)
	       */
	      else
		{
		  strcpy (tokenizer, struct_Ptr->loc_list[k]);

		  if (tokenizer[0] != "\""[0])
		    {
		      pch = strtok_r (tokenizer, "\"", &dump_calc);
		      strcpy (buffer, pch);
		      pch = strtok_r (NULL, "\"", &dump_calc);
		    }
		  else
		    {
		      strcpy (buffer, "");
		      pch = strtok_r(tokenizer, "\"", &dump_calc);
		      pch++;
		    }

		  if (struct_Ptr->vars_types[struct_Ptr->loc_vars_ranks
					     [k][0]]
		      == DISCRETE)
		    cache[struct_Ptr->loc_vars_ranks[k][0]].num_value
		      = atof(cache
			     [struct_Ptr->loc_vars_ranks[k][0]]
			     .string_value);

		  sprintf (buffer, "%s%f", buffer,
			   cache[struct_Ptr->loc_vars_ranks[k][0]]
			   .num_value);

		  for (p = 1; p < struct_Ptr->loc_vars_count[k]; ++p)
		    {
		      pch = strtok_r (NULL, "\"", &dump_calc);
		      strcat (buffer, pch);
		      pch = strtok_r (NULL, "\"", &dump_calc);

		      sprintf (buffer, "%s%f", buffer,
			       cache[struct_Ptr->loc_vars_ranks[k][p]]
			       .num_value);
		    }

		  pch = strtok_r (NULL, "\"", &dump_calc);

		  if (pch != NULL)
		    strcat(buffer, pch);

		  /* appel d'un équivalent C du 'eval' pythonesqe */
		  return_check = evaluator.evaluate (buffer, &num_value);

		  if (return_check == -2)
		    {
		      printf("Erreur de champ (range) dans les calculs \
locaux: %s, calcul: %s, scénario: %s, iteration %d\n",
			     struct_Ptr->loc_list[k], buffer,
			     struct_Ptr->name, i);
		      exit(1);
		    }

		  cache [struct_Ptr->vars_count + struct_Ptr->c_bool_count
			 + k].num_value = num_value;
		  struct_Ptr->loc_results [offset_loc + k] += num_value ;
		}
	    }
	}
      gzclose (p_file);

      /* Pour pouvoir afficher une progression */
      struct_Ptr->progress[struct_Ptr->progress_id]++;
    }
  free (cache);

  /* On utilise un pointeur de type void (seul retour possible d'une
   * fonction passée à un thread) pour contenir et retourner un int.
   * Conversion intermédiare pour éviter un avertissement du compilateur.
   */
  return (void *) CONVERSION struct_Ptr->num_scen ;
}

/*
 * Chargé de mettre à jour et d'afficher la progression du parsing à
 * intervalle régulier.
 */
void * display_progress (void *ptr)
{
  /* Transformer un pointeur void en pointeur de struct */
  pBar_args *current_Prog = (pBar_args*) ptr;

  /* Point d'entrée */
  fpos_t orig;
  fgetpos (current_Prog->progress_file, &orig);

  fprintf (current_Prog->progress_file, "0.00\n");
  struct timespec wait_time = {SLEEP_TIME, 0};

  float percent_progress = 0;
  while (percent_progress != 100)
    {
      percent_progress = 0;

      /* Aucune pénalité en performance à simplement faire du polling */
      nanosleep (&wait_time, NULL);
      for (int i = 0; i < current_Prog->threads_count; i++)
	percent_progress += current_Prog->progress [i];

      percent_progress *= 100. / current_Prog->total_iters;

      /* Revenir au début du fichier */
      fsetpos (current_Prog->progress_file, &orig);

      fprintf (current_Prog->progress_file, "%.2f\n", percent_progress);
    }
  return NULL;
}

/*
 * Affiche les résultats du parsing et fait les calculs demandés.
 */
void print_results (print_func_args *args)
{
  /* Les différents offsets dûs au scénario en cours */
  int offset_acc  = args->num_scen * args->iters_count
    * args->acc_vars_count;
  int offset_bo   = args->num_scen * args->iters_count
    * args->bool_vars_count;
  int offset_dis  = args->num_scen * args->iters_count
    * args->discrete_vars_count;
  int offset_loc  = args->num_scen * args->iters_count
    * args->total_loc_count;
  int offset_c_bo = args->num_scen * args->iters_count * args->c_bool_count;
  int icr_off     = args->num_scen * args->ICR_vars_count;

  puts("---------------------------------------");
  printf("Nom du scénario: %s\n", args->name);
  printf("Population: %d\n", args->pop);
  printf("Nombre d'itérations: %d\n", args->iters_count);
  puts("---------------------------------------\n");
  puts("Variables standards:");
  puts("-----------------------------\n");
  puts("Désignation,Type,Valeur,Relatif (% ou par individu),Ecart type,\
IC (±)");

  double mean, std, sum;
  int i, v, p;

  set <string> keys;
  set <string>::iterator key_it;
  unordered_map <string, unsigned int>::iterator it;

  int acc_vars_rank  = 0;
  int bool_vars_rank = 0;
  int dis_vars_rank  = 0;

  /* Moyenne + écart type + sauvegarder valeur si dans les ICRs */
  for (i = 0; i < args->vars_count; ++i)
    {
      /* Les variables standards */
      switch (args->vars_types[i])
	{
	case BOOLEAN:
	  /* Pas de sous-routines parce que trop de paramètres à passer */
	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += args->bool_results [offset_bo +
					 (args->bool_vars_count * v)
					 + bool_vars_rank] ;
	    }
	  mean = sum / args->iters_count;

	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += pow (args->bool_results [offset_bo +
					      (args->bool_vars_count * v) +
					      bool_vars_rank] - mean, 2);
	    }
	  std = sqrt (sum / (args->iters_count-1));

	  if (! args->no_show [i] )
	    {
	      printf("%s,Booléenne,%.8G,%.8G %%,%.8G,%.8G\n",
		     args->vars_list[i], mean, (mean / args->pop) * 100,
		     std,  get_CI (std, args->iters_count));
	    }

	  ++bool_vars_rank;
	  break;

	case ACCUMUL:
	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += args->acc_results [offset_acc + (args->acc_vars_count
						      * v) + acc_vars_rank];
	    }
	  mean = sum / args->iters_count;

	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += pow (args->acc_results [offset_acc +
					     (args->acc_vars_count * v) +
					     acc_vars_rank] - mean, 2);
	    }
	  std = sqrt (sum / (args->iters_count-1));

	  if (! args->no_show [i] )
	    {
	      printf("%s,Accumulatrice,%.8G,%.8G,%.8G,%.8G\n",
		     args->vars_list[i], mean, mean / args->pop, std,
		     get_CI (std, args->iters_count));
	    }

	  ++acc_vars_rank;
	  break;

	case DISCRETE:
	  /* On insère d'abord dans le « set » les clés existantes */
	  for (v = 0; v < args->iters_count; ++v)
	    for (it = args->discrete_results[offset_dis +
					     (args->discrete_vars_count * v)
					     + dis_vars_rank].begin();
		 it != args->discrete_results[offset_dis +
					      (args->discrete_vars_count
					       * v) + dis_vars_rank].end();
		 ++it)
	      keys.insert(it->first);

	  if (! args->no_show [i])
	    printf("%s, Discrète\n", args->vars_list[i]);

	  /* Toutes les valeurs pour les clés existantes */
	  for (key_it = keys.begin(); key_it != keys.end(); ++key_it)
	    {
	      sum = 0;
	      for (v = 0; v < args->iters_count; ++v)
		{
		  sum += args->discrete_results[offset_dis +
						(args->discrete_vars_count
						 * v) + dis_vars_rank]
		    [*key_it];
		}
	      mean = sum / args->iters_count;

	      sum = 0;
	      for (v = 0; v< args->iters_count; ++v)
		{
		  sum += pow (args->discrete_results
			      [offset_dis + (args->discrete_vars_count
					     * v) + dis_vars_rank]
			      [*key_it] - mean, 2);
		}
	      std = sqrt (sum / (args->iters_count-1));

	      if (! args->no_show [i] )
		{
		  printf(",%s,%.8G, %.8G %%, %.8G, %.8G\n", key_it->c_str(),
			 mean, (mean / args->pop) * 100, std,
			 get_CI (std, args->iters_count));
		}
	    }

	  ++dis_vars_rank;
	  keys.clear(); /* Remettre à zéro pour itération suivante */
	  break;
	}

      if (! args->ICR_vars_count)
	continue;

      if (args->vars_types[i] < DISCRETE)
	{
	  /* Sauvegarder certaines valeurs pour les ICR */
	  if (args->cmp_type < DISCRETE && args->cmp_rank == i)
	    {
	      args->cmp_means [args->num_scen] = mean;
	      args->cmp_stds  [args->num_scen] = std;
	    }
	  else
	    {
	      for (v = 0; v < args->ICR_vars_count; ++v)
		{
		  if (args->ICR_vars_types [v]
		      < DISCRETE && args->ICR_vars_ranks [v] == i)
		    {
		      args->ICR_vars_means [icr_off + v] = mean;
		      args->ICR_vars_stds  [icr_off + v] = std;
		      break;
		    }
		}
	    }
	}
    }

  /* Expressions booléennes */
  if (args->c_bool_count)
    {
      puts("\n\nExpressions booléennes:");
      puts("-----------------------------\n");
      puts("Désignation,Valeur,Relatif,Ecart type,IC (±)");

      for (i = 0; i < args->c_bool_count; i++)
	{
	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += args->c_bool_results [offset_c_bo +
					   (args->c_bool_count * v) + i];
	    }
	  mean = sum / args->iters_count;

	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += pow (args->c_bool_results [offset_c_bo +
						(args->c_bool_count * v)
						+ i] - mean, 2);
	    }
	  std = sqrt (sum / (args->iters_count-1));

	  if (! args->no_show [args->vars_count + i] )
	    {
	      printf("%s,%.8G,%.8G %%,%.8G,%.8G\n", args->c_bool_labels[i],
		     mean, (mean / args->pop) * 100, std,
		     get_CI (std, args->iters_count));
	    }

	  if (! args->ICR_vars_count)
	    continue;

	  /* Sauvegarder certaines valeurs pour les ICR */
	  if (args->cmp_type == CUSTOM_BOOLEAN && args->cmp_rank == i)
	    {
	      args->cmp_means [args->num_scen] = mean;
	      args->cmp_stds  [args->num_scen] = std;
	    }
	  else
	    {
	      for (v = 0; v < args->ICR_vars_count; ++v)
		{
		  if (args->ICR_vars_types [v] == CUSTOM_BOOLEAN
		      && args->ICR_vars_ranks [v] == i)
		    {
		      args->ICR_vars_means [icr_off + v] = mean;
		      args->ICR_vars_stds  [icr_off + v] = std;
		      break;
		    }
		}
	    }
	}
    }

  /* Calculs locaux */
  if (args->total_loc_count)
    {
      char *cond_var;
      puts("\n\nCalculs locaux:");
      puts("-----------------------------\n");
      puts("Désignation,Calcul effectué,Condition,Valeur,Ecart type,IC (±)");

      for (i = 0; i < args->total_loc_count; ++i)
	{
	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += args->loc_results [offset_loc + (args->total_loc_count
						      * v) + i] ;
	    }
	  mean = sum / args->iters_count;

	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += pow (args->loc_results [offset_loc
					     + (args->total_loc_count
						* v) + i]
			  - mean, 2) ;
	    }
	  std = sqrt (sum / (args->iters_count-1));

	  if (! args->no_show [args->vars_count + args->c_bool_count + i] )
	    {
	      if (args->loc_labels[i] != NULL)
		printf("%s,", args->loc_labels[i]);
	      else
		printf(",");

	      /* Vérifier s'il s'agit d'un calcul avec condition */
	      if ( i >= args->loc_count )
		{
		  cond_var = args->cond_vars_rank[i - args->loc_count]
		    >= args->vars_count ?
		    args->c_bool_labels[args->cond_vars_rank
					[i- args->loc_count]
					- args->loc_count
					- args->vars_count]
		    : args->vars_list[args->cond_vars_rank
				     [i - args->loc_count]];

		  printf("%s,%s", args->loc_list[i], cond_var);
		}
	      else
		printf("%s,", args->loc_list[i]);

	      printf(",%.8G,%.8G,%.8G\n", mean, std,
		     get_CI (std, args->iters_count));
	    }

	  if (! args->ICR_vars_count)
	    continue;

	  /* Sauvegarder certaines valeurs pour les ICR */
	  if (args->cmp_type == LOC_CALC && args->cmp_rank == i)
	    {
	      args->cmp_means [args->num_scen] = mean;
	      args->cmp_stds  [args->num_scen] = std;
	    }
	  else
	    {
	      for (v = 0; v < args->ICR_vars_count; ++v)
		{
		  if (args->ICR_vars_types [v] == LOC_CALC
		      && args->ICR_vars_ranks [v] == i)
		    {
		      args->ICR_vars_means [icr_off + v] = mean;
		      args->ICR_vars_stds  [icr_off + v] = std;
		      break;
		    }
		}
	    }
	}
    }

  /* Calculs globaux */
  if (args->calcs_count)
    {
      char   tokenizer [BUFFER_SIZE]; /* pour recevoir une copie qui sera
                                       * mise en morceaux par strtok_r */
      char   buffer    [BUFFER_SIZE];
      char   *parsing, *dump;
      int    nb_equ_vars;

      double *storing = (double*) malloc
	(sizeof(double) * args->iters_count * args->calcs_count);

      int    return_check;
      int    glob_offs;
      int    premature_exit;

      eval   evaluator;

      puts("\n\nCalculs globaux:");
      puts("-----------------------------\n");
      puts("Désignation,Calcul effectué,Valeur,Ecart type,IC (±)");

      /* Même procédure que calculs locaux. Un peu plus complexe car
       * + de types possibles. (géré par des énoncés "switch") */
      for (i = 0; i < args->calcs_count; ++i)
	{
	  glob_offs = i * args->iters_count;
	  premature_exit = 0;
	  sum = 0;
	  parsing = strchr (args->calcs_list[i], "\""[0]);

	  if (parsing == NULL)
	    {
	      printf("Calcul invalide! Il faut utiliser au moins une \
variable: %s\n",
		     args->calcs_list[i]);
	      continue;
	    }

	  nb_equ_vars = 0;

	  while (parsing != NULL)
	    {
	      parsing = strchr (parsing+1, "\""[0]);
	      ++nb_equ_vars;
	    }
	  nb_equ_vars >>= 1;

	  for (v = 0; v < args->iters_count; ++v)
	    {
	      strcpy (tokenizer, args->calcs_list[i]);

	      if (tokenizer[0] != "\""[0])
		{
		  parsing = strtok_r (tokenizer, "\"", &dump);
		  strcpy (buffer, parsing);
		  parsing = strtok_r (NULL, "\"", &dump);
		}
	      else
		{
		  strcpy(buffer, "");
		  parsing = strtok_r(tokenizer, "\"", &dump);
		  parsing++;
		}

	      switch (args->calcs_vars_types[i][0])
		{
		case BOOLEAN:
		  sprintf (buffer, "%s%u", buffer, args->bool_results
			   [offset_bo + (args->bool_vars_count * v) +
			    args->calcs_relative_ranks[i][0]]);
		  break;

		case ACCUMUL:
		  sprintf (buffer, "%s%f", buffer, args->acc_results
			   [offset_acc + (args->acc_vars_count * v)
			    + args->calcs_relative_ranks[i][0]]);
		  break;

		case CUSTOM_BOOLEAN:
		  sprintf (buffer, "%s%u", buffer, args->c_bool_results
			   [offset_c_bo + (args->c_bool_count * v)
			    + args->calcs_relative_ranks[i][0]]);
		  break;

		case LOC_CALC:
		  sprintf (buffer, "%s%f", buffer, args->loc_results
			   [offset_loc + (args->total_loc_count * v)
			    + args->calcs_relative_ranks[i][0]]);
		  break;

		case GLOB_CALC:
		  sprintf (buffer, "%s%f", buffer, storing
			   [(args->calcs_relative_ranks[i][0]
			     * args->iters_count) + v]);
		  break;
		}

	      for (p = 1; p < nb_equ_vars; ++p)
		{
		  parsing = strtok_r (NULL, "\"", &dump);
		  strcat (buffer, parsing);
		  parsing = strtok_r (NULL, "\"", &dump);

		  switch (args->calcs_vars_types[i][p])
		    {
		    case BOOLEAN:
		      sprintf (buffer, "%s%u", buffer, args->bool_results
			       [offset_bo + (args->bool_vars_count * v)
				+ args->calcs_relative_ranks[i][p]]);
		      break;

		    case ACCUMUL:
		      sprintf (buffer, "%s%f", buffer, args->acc_results
			       [offset_acc + (args->acc_vars_count * v)
				+ args->calcs_relative_ranks[i][p]]);
		      break;

		    case CUSTOM_BOOLEAN:
		      sprintf (buffer, "%s%u", buffer, args->c_bool_results
			       [offset_c_bo + (args->c_bool_count * v)
				+ args->calcs_relative_ranks[i][p]]);
		      break;

		    case LOC_CALC:
		      sprintf (buffer, "%s%f", buffer, args->loc_results
			       [offset_loc + (args->total_loc_count * v)
				+ args->calcs_relative_ranks[i][p]]);
		      break;

		    case GLOB_CALC:
		      sprintf (buffer, "%s%f", buffer, storing
			       [(args->calcs_relative_ranks[i][p]
				 * args->iters_count) + v]);
		      break;
		    }
		}

	      parsing = strtok_r (NULL, "\"", &dump);

	      if (parsing != NULL)
		strcat(buffer, parsing);

	      /* appel d'un équivalent C du 'eval' pythonesqe */
	      return_check = evaluator.evaluate (buffer, storing
						 + glob_offs + v);

	      if (return_check == -2)
		{
		  printf("Erreur de champ (« range ») dans les calculs \
globaux: %s, calcul: %s, scénario: %s, iteration: %d\n",
			 args->calcs_list[i], buffer, args->name, v);
		  premature_exit = 1;

		  if (args->calcs_labels[i] != NULL)
		    printf("Avertissement: Une étiquette a été détectée: \
%s. Prendre garde que toute expression/calcul utilisant cette \
variable sera évidemment faussé(e).\n",
			   args->calcs_labels[i]);
		}

	      sum += storing [glob_offs + v];
	    }

	  mean = sum / args->iters_count;

	  sum = 0;
	  for (v = 0; v < args->iters_count; ++v)
	    {
	      sum += pow (storing[glob_offs + v] - mean, 2);
	    }
	  std = sqrt(sum / (args->iters_count-1));

	  if (! args->no_show [args->vars_count + args->c_bool_count
			       + args->total_loc_count + i]
	      && ! premature_exit)
	    {
	      if (args->calcs_labels[i] != NULL)
		printf("%s,", args->calcs_labels[i]);
	      else
		printf(",");

	      printf("%s,%.8G,%.8G,%.8G\n", args->calcs_list[i], mean, std,
		     get_CI (std, args->iters_count));
	    }

	  if (! args->ICR_vars_count)
	    continue;

	  /* Sauvegarder certaines valeurs pour les ICR */
	  if (args->cmp_type == GLOB_CALC && args->cmp_rank == i)
	    {
	      args->cmp_means [args->num_scen] = mean;
	      args->cmp_stds  [args->num_scen] = std;
	    }
	  else
	    {
	      for (v = 0; v < args->ICR_vars_count; ++v)
		{
		  if (args->ICR_vars_types [v] == GLOB_CALC
		      && args->ICR_vars_ranks [v] == i)
		    {
		      args->ICR_vars_means [icr_off + v] = mean;
		      args->ICR_vars_stds  [icr_off + v] = std;
		      break;
		    }
		}
	    }
	}
      free (storing);
    }
  printf("\n\n");
  return;
}

/*
 * Callback de qsort.
 */
int sort_func (const void *elem1, const void *elem2)
{
  return ( *(double*)elem1 - *(double*)elem2 ) > 0 ? 1 : -1 ;
}
