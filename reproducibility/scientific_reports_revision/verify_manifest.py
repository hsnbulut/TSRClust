"""Verify published file identities without executing R or submitting jobs."""
from pathlib import Path
import hashlib,csv
D=Path(__file__).resolve().parent
rows=list(csv.DictReader((D/'SHA256SUMS.csv').open()))
for row in rows:
 p=D/row['path'];assert p.is_file() and p.stat().st_size==int(row['bytes'])
 assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
print('PUBLIC_MANIFEST_OK',len(rows),'files')
