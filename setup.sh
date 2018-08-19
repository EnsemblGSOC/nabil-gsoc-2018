sudo -H pip3 install -r requirements.txt
cd weighted_alignment
g++ weighted_needleman_wunsch.cpp -o exe_file
cd ..
echo "----------- Setup Completed -----------"