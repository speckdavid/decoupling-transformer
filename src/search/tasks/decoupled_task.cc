#include "decoupled_task.h"

#include "factoring.h"

#include "../task_utils/task_dump.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace extra_tasks {
DecoupledTask::DecoupledTask(const shared_ptr<AbstractTask> &parent)
    : DelegatingTask(parent) {
        Factoring factoring(parent);

    ofstream output_file; // Creating an output file stream object

    // Open a file named "example.txt" for writing
    output_file.open("decoupled_task.sas");
    task_dump::dump_as_SAS(*parent, output_file);
    output_file.close();
}
}
