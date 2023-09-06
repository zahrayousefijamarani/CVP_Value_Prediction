
for f in ../public_kit/compute_int/*
do
./cvp -v -t 1 -F 16,0,0,0,0 -f 5 -M 0 -A 0 -w 256 -D 15,4,64,2,20,8,64,12,23,16,128,60,150 "$f" >> "$f.txt" 
echo "done"
done

for f in ../public_kit/compute_fp/*
do
./cvp -v -t 1 -F 16,0,0,0,0 -f 5 -M 0 -A 0 -w 256 -D 15,4,64,2,20,8,64,12,23,16,128,60,150 "$f" >> "$f.txt"
done

for f in ../public_kit/srv/*
do
./cvp -v -t 1 -F 16,0,0,0,0 -f 5 -M 0 -A 0 -w 256 -D 15,4,64,2,20,8,64,12,23,16,128,60,150 "$f" >> "$f.txt"
done


