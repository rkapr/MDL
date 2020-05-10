mkdir edited_4;
for file in random_network4_*csv; do awk '$0=$0"];"' $file > edited_4/temp.csv;sed "s/,/=[/" edited_4/temp.csv > edited_4/$file; mv edited_4/$file edited_4/${file%.csv}.m ; done
