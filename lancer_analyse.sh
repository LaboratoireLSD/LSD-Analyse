#!/bin/bash

if [ $# -lt 1 ]
then
    echo -e "\033[1mUsage:\033[0m `basename $0` répertoire-cible \
\033[2m[fichier-configuration]\033[0m"
    echo -e "\033[1mAide:\033[0m `basename $0` -? / -h / --help"
    exit 0
else
    if [ $1 = "-?" -o $1 = "-h" -o $1 = "--help" ]
    then
	man `basename $0` || exit 1
	exit 0
    fi
fi

main_dir=$1/
dir_results=$1/Results/
dir_analyse=$1/Analyse/

if [ ! -d "${main_dir}" ]
then
    echo "Le répertoire spécifié est inexistant: ${main_dir}"
    exit 1
fi

if [ $# -gt 1 ]
then
    if [ -f "$2" ]
    then
        conf=$2
    else
        echo "Fichier de configuration inexistant: $2"
        exit 1
    fi
else
    conf=${dir_analyse}default.conf

    if [ ! -f "${conf}" ]
    then
        echo -e "; Ceci est un template de fichier de configuration\n\n\
[proportions]\n\n\n[calculs (local)]\n\n\n[expressions booleennes]\n\n\n\
[calculs (conditionnel)]\n\n\n[calculs (global)]\n\n\n[ICER]\n\n\n\
[ne pas afficher]\n" > ${dir_analyse}default.conf
    fi
fi

recode dos..lat1 -q ${conf}

uni=0
for sub in "${main_dir}univariate"*
do
    if [ -d "${sub}" ]
    then
	$0 ${sub} ${conf}
	uni=1
    fi
done

if [ ! -d "${dir_results}" ]
then
    if (( uni )); then
	exit 0
    else
	echo "Le répertoire de résulats est inexistant: ${dir_results}"
	exit 1
    fi
fi

mkdir -p ${dir_analyse} || exit 1

resultats=$(basename ${conf})
resultats=${resultats%"."*}.txt

rm -f ${dir_analyse}${resultats} || exit 1

if (( $(find ${main_dir} -type f -name 'config-*.tar.gz' | wc -l) ))
then
    cd ${main_dir}
    tarFile=$(ls -t config-*.tar.gz | head -1)

    if [[ $(basename $(pwd)) = "univariate_"* ]]
    then
	nomUnivar=$(basename $(pwd))
	nomUnivar=${nomUnivar:11:-3}
	valUnivar=$(tar -zxOf ${tarFile} ./parameters_0.xml | pcregrep -M \
"${nomUnivar}.*\n.*value=\".*\"" | grep -oE "[0-9]+\.?[0-9]*")

	echo -e "-!- Analyse univariée détectée -!-\nValeur actuelle de la \
variable '${nomUnivar}': ${valUnivar}\n" | tee Analyse/${resultats} || exit 1
    else
	tar -zxf ${tarFile} ./parameters.xml
	mv -f parameters.xml Analyse
    fi
    cd - > /dev/null
fi

arrayScenarios=()
for scenario in "${dir_results}"*
do
    if [ -d "${scenario}" ]
    then
	arrayScenarios+=(${scenario:${#dir_results}})
    fi
done

nbSims=$(find ${dir_results}${arrayScenarios[0]} -type f -name *_Output.gz \
| wc -l)

nbThreads=$(grep -c processor /proc/cpuinfo)

{
    Analyse ${main_dir} ${conf} ${nbSims} ${nbThreads} \
${arrayScenarios[@]} || echo -e "Le programme C++ a été interrompu prématurément."
}   | tee -a ${dir_analyse}${resultats} || exit 1

sed -i -e 's/$/\r/' ${dir_analyse}${resultats} ${conf}

if [[ $(tail -1 ${dir_analyse}${resultats}) =~ "Le programme C++ a été interrompu prématurément." ]]
then
    exit 1
else
    exit 0
fi

