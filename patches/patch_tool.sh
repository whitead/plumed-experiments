#!/bin/bash
# This script is not executed directly, but just sourced by all the code-specific scripts.
# It is meant to be as general as possible, to simply apply the directives set in the code-specific scripts,
# and to take care of the patch. In particular, it keeps track of applied patches, to avoid re-patching
# and other errors. It is rather verbose, so it should be easy to follow the way it works simply reading
# the messages.
# In short, it can do three tasks:
# -patch: It applies the patch to an out-of-the-box code. Changed files are backed up adding the extension .preplumed
# -revert: It reverts the code. To do that, it uses the backup files (.preplumed)
# -save: (for developers) It looks for .preplumed files and save the difference between them and the modified ones
#        in a patch script in the patches/diff directory. This script will be then used to -patch.
#
# The idea behind the save task is that one directly works on the original files, then save the patch.
# If you have to modify extra files, just copy them (cp pippo.c pippo.c.preplumed) and the -save option will
# automatically take care of them at the next save. Also, if you are working on the CVS version, don't forget
# to commit di files in the patches/diff directory

PLUMED_VERSION="1.3"
DIFF_SCRIPT="$plumedir/patches/diff/$CODE"

# everything is performed in the destination directory
cd "$destination"

if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$NAME  (-patch) (-revert)   "
 echo " -patch  : apply PLUMED patch " 
 echo " -revert : revert code to original "
 echo " (for developers:"
 echo "    -save : saves the changes done to patched files "
 echo " )"
 exit
fi

case "$1" in
(-patch)
  echo "* I will try to patch PLUMED version $PLUMED_VERSION ..."
  if test -e "PLUMED_PATCH" ; then
    echo "-- File $destination/PLUMED_PATCH exists"
    echo "-- I guess you have already patched PLUMED $(tail -1 PLUMED_PATCH)"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
  fi
  echo "#Please do not remove or modify this file"                    >  PLUMED_PATCH
  echo "#It is keeps track of patched versions of the PLUMED package" >> PLUMED_PATCH
  echo "$PLUMED_VERSION"                                              >> PLUMED_PATCH
  
  echo "-- Executing pre script"

  command -v patch &>/dev/null || { echo "I require patch command but it's not installed. Aborting." >&2; exit 1;  } 

#------------------- first, check if GNU patch works
cat > test_patch1 << \EOF
alfa
beta
EOF

cat > test_patch2 << \EOF
alfa
gamma
EOF

cat > test_patch3 << \EOF_EOF
patch -c -l -b -F 3 --suffix=.preplumed "./test_patch1" << \EOF
EOF_EOF

diff -c test_patch1 test_patch2 >> test_patch3

echo EOF >> test_patch3

bash test_patch3 &> test_patch4

status=$?
if [ $status -ne 0 ]
then
  echo "patch does not work! Error message:"
  echo "**********"
  cat test_patch4
  echo "**********"
  echo "Please install a recent version of the GNU patch utility and try again."
  exit
fi

rm test_patch1 test_patch2 test_patch3 test_patch4
if [ -e test_patch1.preplumed ]
then
  rm test_patch1.preplumed
fi
#-------------------

  to_do_before_patch

  echo "-- Setting up symlinks"
  for file in $LINKED_FILES ; do
    base="${file##*/}"
    if test -e $WHERE_LINKS/$base ; then
      echo "PATCH ERROR: file $base is already in $WHERE_LINKS"
      exit 1
    fi
    if test -n "$PLUMED_RENAME_RULE" ; then
      base="$($PLUMED_RENAME_RULE $base)"
    fi
    ln -s $file $WHERE_LINKS/$base
  done

  if [ "$RECON_LIBS" != "" ] ; then
    echo "-- Setting up recon symlinks"
    for file in $RECON_LINKS ; do
      base="${file##*/}"
      if test -e $WHERE_LINKS/$base ; then
        echo "PATCH ERROR: file $base is already in $WHERE_LINKS"
        exit 1
      fi
      if test -n "$RECON_RENAME_RULE" ; then
        base="$($RECON_RENAME_RULE $base)"
      fi
      ln -s $file $WHERE_LINKS/$base
    done
  fi 

  if test -e "$DIFF_SCRIPT" ; then
    echo "-- Applying patches"
    bash "$DIFF_SCRIPT"
  fi

  echo "-- Executing post script"
  to_do_after_patch

  echo "- DONE!"
;;
(-revert)
  echo "* I will try to revert PLUMED version $PLUMED_VERSION ..."
  if test ! -e PLUMED_PATCH ; then
    echo "-- File $destination/PLUMED_PATCH is not there"
    echo "-- I guess you never patched, so there is nothing to revert"
    echo "-- See you later"
    echo "* ABORT"
    exit
  fi
  PATCHED_VERSION="$(tail -1 $destination/PLUMED_PATCH)"
  if test "$PLUMED_VERSION" != "$PATCHED_VERSION" ; then
    echo "-- I guess from file $destination/PLUMED_PATCH"
    echo "   that you patched a different version ($PATCHED_VERSION)"
    echo "-- If you wish to revert, you have to use the same script that you used to patch"
    echo "-- See you later"
    echo "* ABORT"
  fi

  echo "-- Executing pre script"

  to_do_before_revert

  echo "-- Restoring .preplumed files"
  PREPLUMED=$(find . -name "*.preplumed")
  if ! test "$PREPLUMED" ; then
    echo "-- I cannot find any .preplumed file"
    echo "* ABORT"
    exit
  fi

  for bckfile in $PREPLUMED ; do
    file="${bckfile%.preplumed}"
    mv "$bckfile" "$file"
  done

  echo "-- Removing symlinks"
  for file in $LINKED_FILES ; do
    base="${file##*/}"
    if test -n "$PLUMED_RENAME_RULE" ; then
      base="$($PLUMED_RENAME_RULE $base)"
    fi
    if test -e $WHERE_LINKS/$base ; then
      rm $WHERE_LINKS/$base 
    else
      echo "PATCH WARNING: file $base is not in $WHERE_LINKS"
    fi
  done

  if [ "$RECON_LIBS" != "" ] ; then
    echo "-- Removing recon symlinks"
    for file in $RECON_LINKS ; do
      base="${file##*/}"
      if test -n "$RECON_RENAME_RULE" ; then
        base="$($RECON_RENAME_RULE $base)"
      fi
      if test -e $WHERE_LINKS/$base ; then
        rm $WHERE_LINKS/$base
      else
        echo "PATCH WARNING: file $base is not in $WHERE_LINKS"
      fi
    done
  fi

  echo "-- Executing post script"
  to_do_after_revert

  rm "$destination/PLUMED_PATCH"

  echo "* DONE!"
;;
(-save)
  echo "* I will try to save the changes done in $destination"
  echo "* THIS IS ONLY FOR DEVELOPERS"
  if test ! -e PLUMED_PATCH ; then
    echo "-- File $destination/PLUMED_PATCH is not there"
    echo "-- To save your changes, you need to work on a patched version"
    echo "* ABORT"
    exit
  fi
  PATCHED_VERSION="$(tail -1 $destination/PLUMED_PATCH)"
  if test "$PLUMED_VERSION" != "$PATCHED_VERSION" ; then
    echo "-- I guess from file $destination/PLUMED_PATCH"
    echo "   that you patched a different version ($PATCHED_VERSION)"
    echo "-- To save your changes, you need to work on the same version"
    echo "* ABORT"
    exit
  fi

  echo "-- First looking for .preplumed files in $destination"
  PREPLUMED=$(find . -name "*.preplumed")
  if ! test "$PREPLUMED" ; then
    echo "-- I cannot find any .preplumed file"
    echo "* ABORT"
    exit
  fi
  echo "-- Found the following:"
  echo "$PREPLUMED"
  echo "-- I build now a diff script in $DIFF_SCRIPT ..."

  test -e "$DIFF_SCRIPT" && rm "$DIFF_SCRIPT"

  for bckfile in $PREPLUMED ; do
    file="${bckfile%.preplumed}"
    if test ! -e "$file" ; then
      echo "-- File $file is not there, thus I cannot take the diff"
      echo "* ABORT"
    else
      echo "patch -c -l -b -F 3 --suffix=.preplumed \"${file}\" << \\EOF_EOF" >> "$DIFF_SCRIPT"
      diff -c "${bckfile}" "$file"                                     >> "$DIFF_SCRIPT"
      echo "EOF_EOF"                                                   >> "$DIFF_SCRIPT"
    fi
  done

  echo "-- ... done"
  echo "-- The best way to check it is to revert, patch again and see if everything is in the proper place"
  echo "-- Bye!"
  
esac


