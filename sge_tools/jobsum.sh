qacct -j $1 > /tmp/$1.acct
grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | awk '{s+=$1} END {print "Total CPU time : " s " sec"}'
grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | awk '{s+=$1} END {print "mean CPU time : " s/NR " sec"}'
echo Min CPU time : $(grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | sort -n | head -n1) sec
echo Max CPU time : $(grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | sort -n | tail -n1) sec
submittime=`grep qsub_time /tmp/$1.acct | head -1 | sed "s/qsub_time    //"`
finishtime=`grep end_time /tmp/$1.acct | tail -1 | sed "s/end_time    //"`
echo Submitted at: $submittime
echo Finished at: $finishtime
submit_epoch=`date -d "$submittime" +%s`
finished_epoch=`date -d "$finishtime" +%s`
echo Submit-to-finish time : $((finished_epoch - submit_epoch)) sec 
echo Max max memory : $(grep maxvmem /tmp/$1.acct | sort -n | tail -n1 | sed "s/maxvmem    //")
