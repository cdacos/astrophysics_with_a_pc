import sys

def start_parameter(text, i):
  if len(sys.argv) > i:
    print('{0}{1}'.format(text, sys.argv[i]))
    return float(sys.argv[i])
  else:
    return float(raw_input(text))
