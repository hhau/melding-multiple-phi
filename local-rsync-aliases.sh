# TODO: get these names dynamically?

alias push-latest-code=\
"
rsync \
  -avr \
  -e ssh \
  --update \
  --exclude rds \
  --exclude plots \
  --exclude logs \
  . \
  aam71@login.hpc.cam.ac.uk:TEMPLATENAME
" 

alias pull-latest-results=\
"
rsync \
  -avr \
  -e ssh \
  --update \
  aam71@login.hpc.cam.ac.uk:TEMPLATENAME/ \
  . 
" 