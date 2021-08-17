for run in $(ls | grep GB21) ; do prefix=$(echo $run | sed 's/_f.*$//') ; ../nanosims-processor.py --roi $run/${run}_ROI.png -i Aunused/$prefix.im -o $run/$prefix; done
