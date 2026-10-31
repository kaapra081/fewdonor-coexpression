import os, sys, subprocess
log = sys.argv[1]; cmd = sys.argv[2:]
subprocess.Popen(cmd, stdout=open(log,'a'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True, close_fds=True)
