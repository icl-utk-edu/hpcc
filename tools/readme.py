import re
import sys

if len(sys.argv) > 1:
  filename = sys.argv[1]
else:
  filename = "README.html"

d=open(filename).read().replace("</style>","tt {color:navy}\nh2,h3,h4 {color:#527bbd;}\nh2 {border-bottom: 2px solid silver;}\n</style>")
open(filename,"w").write(re.sub("^[.]c000.*",".c000{font-family:monospace;color:navy;}", d, flags=re.MULTILINE))
