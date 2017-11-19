#!/usr/bin/python
import shlex, subprocess

cmd = 'git log -n 1 --date=short --format=format:"rev.%ad.%h" HEAD'
args = shlex.split(cmd)
p = subprocess.Popen(args, stdout=subprocess.PIPE)
sha = p.stdout.readline()
p.terminate()

print "#ifndef VERSION"
print '#define VERSION "%s"' % sha
print "#endif"
