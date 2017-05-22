rm -rf dist build
make CONF=Release build
\rm -f ~/bin/fconv
ln -s $PWD/dist/Release/GNU-Linux-x86/fconv ~/bin/
