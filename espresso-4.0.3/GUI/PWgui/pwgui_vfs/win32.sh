#!/bin/sh

ln -sf make.win32 make.sys
make
zip      pwgui-win32.zip pwgui.exe
tar zcvf pwgui-win32.tgz pwgui.exe
