#include "tools/tool_funcs.h"
#include "tools/competition_tools.h"
#include "gui/gui_manager.h"
#include <math.h>

#ifdef _WIN32
#include <windows.h>  // for SetPriorityClass(), SetThreadPriority()
#endif


// ------------------------------------------------------------

void set_low_priority() {
	// Low thread priority
	#ifdef _WIN32
		if (!SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS)) {
			printf("Failed to set priority on Windows\n");
		}
	#else
		int current_nice = nice(10);  // Increase the nice value (lower priority)
		if (current_nice == -1) {
			perror("Failed to set nice value");
		}
	#endif
}

int main() {
	set_low_priority();

	init_available_systems("../Celestial_Systems/");
	
	print_celestial_system(get_available_systems()[0]);
	
	int selection;
    char title[] = "CHOOSE PROGRAM:";
    char options[] = "Exit; GUI; Run from file";
    char question[] = "Program: ";

    do {
        selection = user_selection(title, options, question);

        switch (selection) {
        case 1:
			start_gui("../GUI/GUI.glade");
            break;
        case 2:
			run_competition_calc("../Queue/example.txt", "../Itineraries/test.itins", get_available_systems()[0]);
            break;
        default:
            break;
        }
    } while(selection != 0);
	

	free_all_celestial_systems();
    return 0;
}