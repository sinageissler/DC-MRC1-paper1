    #!/bin/bash
    #PBS -S /bin/bash
    #PBS -N 2n-rec
    #PBS -o 2n-rec.out
    #PBS -e 2n-rec.err
     
    #PBS -q nogpu_40c_hn
    #PBS -l nodes=1:ppn=20
    #PBS -l walltime=12:00:00
    #PBS -A simlab_project
     
    #PBS -m abe
    #PBS -M geissler@ibpc.fr
     
    #PBS -l epilogue=/shared/scripts/ADMIN__epilogue-qsub.example
     
    ### FOR EVERYTHING BELOW, I ADVISE YOU TO MODIFY THE USER-part ONLY ###
    WORKDIR="/"
    NUM_NODES=$(cat $PBS_NODEFILE|uniq|wc -l)
    if [ ! -n "$PBS_O_HOME" ] || [ ! -n "$PBS_JOBID" ]; then
            echo "At least one variable is needed but not defined. Please touch your manager about."
            exit 1
    else
            if [ $NUM_NODES -le 1 ]; then
                    WORKDIR+="scratch/"
                    export WORKDIR+=$(echo $PBS_O_HOME |sed 's#.*/\(home\|workdir\)/\(.*_team\)*.*#\2#g')"/$PBS_JOBID/"
                    mkdir $WORKDIR
                    rsync -ap $PBS_O_WORKDIR/ $WORKDIR/
            else
                    # WORKDIR+="scratch-dfs/"
                    export WORKDIR=$PBS_O_WORKDIR
            fi       
     
            # if you need to check your job output during execution (example: each hour) you can uncomment the following line
            # /shared/scripts/ADMIN__auto-rsync.example 3600 &
    fi
     
    echo "your current dir is: $PBS_O_WORKDIR"
    echo "your workdir is: $WORKDIR"
    echo "number of nodes: $NUM_NODES"
    echo "number of cores: "$(cat $PBS_NODEFILE|wc -l)
    echo "your execution environment: "$(cat $PBS_NODEFILE|uniq|while read line; do printf "%s" "$line "; done)
     
    cd $WORKDIR
     
    # If you're using only one node, it's counterproductive to use IB network for your MPI process communications
    if [ $NUM_NODES -eq 1 ]; then
            export PSM_DEVICES=self,shm
            export OMPI_MCA_mtl=^psm
            export OMPI_MCA_btl=shm,self
    else
    # Since we are using a single IB card per node which can initiate only up to a maximum of 16 PSM contexts 
    # we have to share PSM contexts between processes
    # CIN is here the number of cores in node
            CIN=$(cat /proc/cpuinfo | grep -i processor | wc -l)
            if [ $(($CIN/16)) -ge 2 ]; then
                    PPN=$(grep $HOSTNAME $PBS_NODEFILE|wc -l)
                    if [ $CIN -eq 40 ]; then
                            export PSM_SHAREDCONTEXTS_MAX=$(($PPN/4))
                    elif [ $CIN -eq 32 ]; then
                            export PSM_SHAREDCONTEXTS_MAX=$(($PPN/2))
                    else
                            echo "This computing node is not supported by this script"
                    fi
                    echo "PSM_SHAREDCONTEXTS_MAX defined to $PSM_SHAREDCONTEXTS_MAX"
            else
    	        echo "no PSM_SHAREDCONTEXTS_MAX to define"
            fi
    fi
     
    function get_gpu-ids() {
    	if [ $PBS_NUM_PPN -eq $(cat /proc/cpuinfo | grep -cE "^processor.*:") ]; then
    		echo "0,1" && return
    	fi
     
    	if [ -e /dev/cpuset/torque/$PBS_JOBID/cpus ]; then
    		FILE="/dev/cpuset/torque/$PBS_JOBID/cpus"
    	elif [ -e /dev/cpuset/torque/$PBS_JOBID/cpuset.cpus ]; then
    		FILE="/dev/cpuset/torque/$PBS_JOBID/cpuset.cpus"
    	else
    		FILE=""
    	fi
     
    	if [ -e $FILE ]; then
    		if [ $(cat $FILE | sed -r 's/^([0-9]).*$/\1/') -eq 0 ]; then
    			echo "0" && return
    		else
    			echo "1" && return
    		fi
    	else
    		echo "0,1" && return
    	fi
    }
     
    gpus=$(get_gpu-ids)
     
    ## END-DO
     
     
    ##
    ## USER part
    ##
     
    ## Environment settings (environment module loadings, etc.)
    module load gcc/9.5.0
    module load gromacs/gnu/2023.2
     
    ## your app calls

    gmx grompp -f step4.0_minimization.mdp -c step3_input.gro -r step3_input.gro -p topol.top -o em.tpr  
    gmx mdrun -v -deffnm em

    equi_prefix=step4.%d_equilibration
    mini_prefix=em

    # run loop of 6 equilibration steps
    cnt=1
    cntmax=6

    while [[ ${cnt} -le ${cntmax} ]]
    do
    pcnt=$((${cnt}-1))
    echo $pcnt
    istep=$(printf ${equi_prefix} ${cnt})
    pstep=$(printf ${equi_prefix} ${pcnt})

    if [[ ${cnt} == 1 ]]
    then
    pstep=${mini_prefix}
    fi 
    gmx grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${pstep}.gro -p topol.top -n index.ndx -maxwarn 3
    gmx mdrun -v -deffnm ${istep} 

    cnt=$((${cnt}+1))
    done

    # Production
    gmx grompp -f step5_production.mdp -c step4.6_equilibration.gro -t step4.6_equilibration.cpt -p topol.top -n index.ndx -o md_500ns.tpr
    #gmx mdrun -deffnm md_500ns 
     
    ## To well chain your jobs, with afterok directive to be sure the current job will complete and OK before running the new one
    # qsub -d `/shared/scripts/getWorkdir.sh` -W depend=afterok:$PBS_JOBID <jobscript>
     
    ##
    ## END-USER part
    ##
     
     
    # At the term of your job, you need to get back all produced data synchronizing workdir folder with you starting job folder 
    # and delete the temporary one (workdir)
    # A good practice is to reduce the file list you need to get back with rsync
    if [ $NUM_NODES -le 1 ]; then
    	cd $PBS_O_WORKDIR
    	rsync -ap $WORKDIR/ $PBS_O_WORKDIR/
    	rm -rf $WORKDIR
    fi
    ## END-DO


