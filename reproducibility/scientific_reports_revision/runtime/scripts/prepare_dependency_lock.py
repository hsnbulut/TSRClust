"""Fetch exact CRAN sources and create a content-hashed, portable source lock."""
from pathlib import Path
import urllib.request,urllib.error,hashlib,tarfile,csv,gzip
root=Path(__file__).resolve().parents[1];base=root/'config/dependencies';cache=base/'sources';cache.mkdir(parents=True,exist_ok=True)
# Preserve submitted clustering engines; choose explicit compatible dependency versions.
pins=[('codetools','0.2-20'),('MASS','7.3-65'),('iterators','1.0.14'),('foreach','1.5.2'),('doParallel','1.0.17'),('DEoptimR','1.1-4'),('Rcpp','1.1.1'),('RcppArmadillo','15.2.4-1'),('rlang','1.2.0'),('generics','0.1.4'),('mvtnorm','1.3-7'),('mclust','6.1.2'),('robustbase','0.99-7'),('dbscan','1.2.4'),('tclust','2.2-0'),('otrimle','2.0'),('teigen','2.2.2')]
rows=[]
for package,version in pins:
 filename=f'{package}_{version}.tar.gz';path=cache/filename
 urls=[f'https://cran.r-project.org/src/contrib/{filename}',f'https://cran.r-project.org/src/contrib/Archive/{package}/{filename}']
 selected=None
 for url in urls:
  try:
   # Always fetch for an intentional new lock; cached lock files are verified by installer.
   data=urllib.request.urlopen(url,timeout=40).read();selected=url;break
  except urllib.error.HTTPError as e:
   if e.code!=404:raise
 if selected is None:raise RuntimeError(f'No CRAN source for {filename}')
 path.write_bytes(data)
 with tarfile.open(path) as tar:
  description=tar.extractfile(f'{package}/DESCRIPTION').read().decode()
 fields={};last=None
 for line in description.splitlines():
  if line.startswith((' ','\t')) and last:fields[last]+=' '+line.strip()
  elif ':' in line:last,val=line.split(':',1);fields[last]=val.strip()
 assert fields['Package']==package and fields['Version']==version
 (base/f'{package}.DESCRIPTION').write_text(description)
 rows.append(dict(package=package,version=version,url=selected,source_file=f'sources/{filename}',sha256=hashlib.sha256(data).hexdigest(),md5=hashlib.md5(data).hexdigest(),bytes=len(data),Depends=fields.get('Depends',''),Imports=fields.get('Imports',''),LinkingTo=fields.get('LinkingTo',''),NeedsCompilation=fields.get('NeedsCompilation','no')))
 print(package,version,selected,flush=True)
with (base/'source_lock.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
