# cimeshare

Repo containing the code from cime/src/share. Used only for moving from cime submodule to directory in E3SM.

How it was made: 
git clone --no-tags git@github.com:ESMCI/cime.git  
(used hash b95a28b417b9b27 from May 4, 2021)

git filter-repo --path-rename src/share:share  
git filter-repo --path share  
git remote add cimeshare (this repo)  
git push cimeshare master  


