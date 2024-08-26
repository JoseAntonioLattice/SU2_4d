BIN=bin
DATA=data

if [ -d $BIN ]
then
   echo $BIN "already exists."
else
    mkdir $BIN
fi


if [ -d $DATA ]
then
   echo $DATA "already exists."
else
    mkdir $DATA
fi
