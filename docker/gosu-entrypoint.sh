#!/bin/bash

# Gosu setup adapted from https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
USER_ID=${RUN_USER_ID:-1000}
GROUP_ID=${RUN_GROUP_ID:-$USER_ID}
USER_NAME=user
GROUP_NAME=viral-ngs

echo "Starting with UID : $USER_ID, GID : $GROUP_ID" >&2
groupadd --gid $GROUP_ID $GROUP_NAME
useradd --shell /bin/bash --uid $USER_ID --gid $GROUP_ID -o -c "" -m $USER_NAME

# Chown to deal with potential problems of creating projects under the path define in the deploy script:
# PROJECTS_PATH="$SCRIPTPATH/$CONTAINING_DIR/$PROJECTS_DIR"
export HOME=/home/$USER_NAME
ln -s /user-data $HOME/data
# This is very slow and not really needed unless /opt/viral-ngs is mounted as a volume on the host
#chown -R $USER_NAME:$GROUP_NAME /opt/viral-ngs

# '-R' is not used here because the user running docker may not have chown -R permission to the whole
# data directory
chown $USER_NAME:$GROUP_NAME /user-data

cd $HOME/data

# Chroot does jump to a custom startup directory, so we don't use it.
# More importantly, it allows to escape to root
# chroot --userspec=$USER_NAME / "/bin/bash"
# Therefore gosu is a good choice
exec /usr/local/bin/gosu $USER_NAME "$@"
