#include <ncurses.h>
#include <locale.h>
#include <wchar.h>
#include <menu.h>

void print_title(){
    setlocale(LC_ALL, "");
    wchar_t* title = L" █████╗ ███████╗████████╗███████╗██████╗ \n██╔══██╗██╔════╝╚══██╔══╝██╔════╝██╔══██╗\n███████║███████╗   ██║   █████╗  ██████╔╝\n██╔══██║╚════██║   ██║   ██╔══╝  ██╔══██╗\n██║  ██║███████║   ██║   ███████╗██║  ██║\n╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝";

    mvaddwstr(1, 1, title);
}

void make_menu(){
    const char* opts[6] = {"Simulation", "Bodies", "Presets", "Advanced", "Fire!", nullptr};
    ITEM* items[6];
    for (int i = 0; i < 5; ++i)
        items[i] = new_item(opts[i], "");

    items[6] = nullptr;

    MENU* menu = new_menu(items);
    post_menu(menu);
    refresh();

    int ch;
    while ((ch = getch()) != '\n') { // till it presses enter
    switch(ch){
        case KEY_DOWN:
            menu_driver(menu, REQ_DOWN_ITEM);
            break;
        case KEY_UP:
            menu_driver(menu, REQ_UP_ITEM);
            break;
    }
    refresh();
    }

    unpost_menu(menu);
    free_menu(menu);
    //for (auto i : items)
    //    free_item(i);
}

int main() {
    initscr();
    start_color();                   
    print_title();
    make_menu();

    init_pair(1, COLOR_RED, COLOR_BLACK);  
    init_pair(2, COLOR_GREEN, COLOR_BLACK);

    refresh();
    getch();
    endwin();
    return 0;
}

