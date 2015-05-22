
import os

Tuple = (
("stream.c", (
  ("HPCC_Stream(", "main"),
  ("checkSTREAMresults(FILE", "checkres"),
  ("checktick() {", "checktick"),
  ("void tuned_STREAM_Copy(", "copy"),
  ("tuned_STREAM_Scale(double", "scale"),
  ("void tuned_STREAM_Add(", "add"),
  ("tuned_STREAM_Triad(double", "triad"),
  ("void computeSTREAMerrors(", "checkerr"),
), "hpcc"),

("stream_mpi.c", (
("main()" , "main"),
("^checktick()", "checktick"),
("computeSTREAMerrors(STREAM", "checkerr"),
("checkSTREAMresults (STREAM", "checkres"),
("void tuned_STREAM_Copy(", "copy"),
("tuned_STREAM_Scale(STREAM", "scale"),
("void tuned_STREAM_Add(", "add"),
("tuned_STREAM_Triad(STREAM", "triad"),
), "tstrm"),
)


def swap_fd(fd, fname, prefix):
  fd.close()

  if not os.path.exists(prefix):
    os.mkdir(prefix)

  name = os.path.join(prefix, fname +".c")
  if fname.startswith("/dev"):
    name = fname
  fd = open(name, "w")
  return fd

        
for tup in Tuple:
  fd = open("/dev/null", "w")

  prefix = tup[2]
  
  for line in open(tup[0]):
    for m in tup[1]:
      if m[0].startswith("^"):
        if line.startswith(m[0][1:]):
          fd = swap_fd(fd, m[1], prefix)

      elif line.find(m[0]) != -1:
        fd = swap_fd(fd, m[1], prefix)

        break

    fd.write(line)

  fd.close()

Replacements = (
("STREAM_TYPE", "double"),
("MAX", "Mmax"),
("MIN", "Mmin"),
("ssize_t", "int"),
("abs", "fabs"),
)

for tup in Tuple:
  prefix = tup[2]
  for m in tup[1]:
    name = m[1]
    if name.startswith("/dev"):
      continue
    fname = os.path.join(prefix, name + ".c")
    code = open(fname).read()
    for rt in Replacements:
      code = code.replace(rt[0], rt[1])
    open(fname, "w").write(code)
