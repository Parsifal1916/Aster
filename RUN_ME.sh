echo "[ ? ] Da non eseguire al di fuori di una cartella"
echo "  *   salva FRAME PER FRAME in file png. Questo verrà"
echo "  *   modificato in futuro e succede perchè mi permette"
echo "  *   di salvare un video completa incollando tutti i frame"
echo
echo "[ * ] Eseguo il compiling..."
g++ -o main src/main.cpp -lsfml-graphics -lsfml-window -lsfml-system -Wno-narrowing -O3 -pthread -std=c++17
echo "[ * ] fatto! eseguibile: ./main"