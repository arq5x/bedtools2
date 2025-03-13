/*****************************************************************************
windowMaker.h

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "NewGenomeFile.h"
#include "bedFile.h"

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class WindowMaker {

public:
    enum INPUT_FILE_TYPE {
        GENOME_FILE,
        BED_FILE
    };
    enum WINDOW_METHOD {
        FIXED_WINDOW_SIZE,
        FIXED_WINDOW_COUNT
    };
    enum ID_METHOD {
        ID_NONE,
        ID_WINDOW_NUMBER,
        ID_SOURCE_ID,
        ID_SOURCE_ID_WINDOW_NUMBER
    };

    // constructor
    WindowMaker(string &fileName, ID_METHOD id_method, INPUT_FILE_TYPE input_file_type, uint32_t count, bool reverse);
    WindowMaker(string &fileName, ID_METHOD id_method, INPUT_FILE_TYPE input_file_type, uint32_t size, uint32_t step, uint32_t initialSize, uint32_t initialStep, bool reverse);

    // destructor
    ~WindowMaker(void);

    void MakeWindowsFromGenome(const string& genomeFileName);
    void MakeWindowsFromBED(string& bedFileName);

private:
    uint32_t _size;
    uint32_t _step;
    uint32_t _count;
    bool _reverse; // should window numbering be reversed?
    WINDOW_METHOD _window_method;
    ID_METHOD _id_method;
    uint32_t _initialSize;
    uint32_t _initialStep;

    void MakeBEDWindow(const BED& interval);

    void MakeFixedSizeWindow(const BED& interval);
    void MakeFixedCountWindow(const BED& interval);
    uint32_t CalculateWindows(const BED& interval);

    string GenerateID(const BED& interval, uint32_t window_index, uint32_t num_windows, bool _reverse) const;
};
