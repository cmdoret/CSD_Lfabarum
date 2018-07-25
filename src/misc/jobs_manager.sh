# This script contains general utilities to manage LSF jobs.
# Cyril Matthey-Doret
# 10.10.2017



function bmonitor {
  # This function hangs the script when the number of queued bjobs
  # containing a given pattern in their name reaches a limit.
  # First argument is the name pattern of jobs to monitor
  # Second argument is the maximum of joobs that can be queued at a time.
  while [ $(bjobs -w | grep "$1" | wc -l) -gt $2 ]
  do
    sleep 1;
  done
}

function prettyload {
  # This function displays a simple animation and shows real-time
  # progress of a task.

  # Spinner animation
  case $(expr $1 % 4) in
    0) sym="|";;
    1) sym="/";;
    2) sym="-";;
    3) sym="\\";;
  esac
  # Computing fill of the progress bar
  progress=$((($1*100)/$2))
  fill=$(printf "%$(($progress/4))s")
  empty=$(printf "%$((25-$progress/4))s")
  printf "   $sym [${fill// /#}${empty// /-}] ${progress}%%"
  # Printing progress
  echo -ne "  $1/$2 \r"
}
