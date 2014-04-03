# grab the first command line argument and use it to call qacct for the job we want to look at
# this takes a few seconds so let's store the output so we can process it multiple times
qacct -j $1 > /tmp/$1.acct

# search for lines that show the cpu time
#                       add up the numbers and print the total as a message
grep cpu /tmp/$1.acct | awk '{s+=$2} END {print "Total CPU time : " s " sec"}'

#                       add up the numbers and divide by the number of lines (NR)
grep cpu /tmp/$1.acct | awk '{s+=$2} END {print "mean CPU time : " s/NR " sec"}'


# ideally we want the max and min CPU time to be similar for efficient cluster use
#                   $() evaluates the commands inside it and outputs the result as text
                                                                           # sort numerically and take the first (smallest)
echo Min CPU time : $(grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | sort -n | head -n1) sec
                                                                           # sort numerically and take the last (biggest)
echo Max CPU time : $(grep cpu /tmp/$1.acct | tr -s ' ' | cut -d ' ' -f2 | sort -n | tail -n1) sec

# to get the submit time...
#             ...search for qsub_time...
#                                          ...take the first line...
#                                                    ...delete the first bit
submittime=$(grep qsub_time /tmp/$1.acct | head -1 | sed "s/qsub_time    //")

# for the finish time...
#                                         we need the end time from the last task   
finishtime=$(grep end_time /tmp/$1.acct | tail -1 | sed "s/end_time    //")
echo Submitted at: $submittime
echo Finished at: $finishtime

# get submit time in UNIX seconds
#              date -d lets you specify a date
#                                    +%s means display the output in UNIX epoch seconds
submit_epoch=$(date -d "$submittime" +%s)
finished_epoch=$(date -d "$finishtime" +%s)

#                               print difference between end and start time in seconds
echo Submit-to-finish time : $((finished_epoch - submit_epoch)) sec 

# we need to keep an eye on the max memory usage
# if this gets too high then we need to increase the number of slots
echo Max max memory : $(grep maxvmem /tmp/$1.acct | sort -n | tail -n1 | sed "s/maxvmem    //")
