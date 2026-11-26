"""Download exactly the existing lock's source archives; never install or rewrite pins."""
from pathlib import Path
import csv,hashlib,urllib.request
D=Path(__file__).resolve().parent/'runtime/config/dependencies'
for row in csv.DictReader((D/'source_lock.csv').open()):
 p=D/row['source_file'];p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():data=p.read_bytes()
 else:
  with urllib.request.urlopen(row['url'],timeout=90) as r:data=r.read()
 assert hashlib.sha256(data).hexdigest()==row['sha256'],row['package']
 if not p.exists():p.write_bytes(data)
 print('SOURCE_PIN_OK',row['package'],row['version'])
