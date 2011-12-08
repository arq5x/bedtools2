/*****************************************************************************
  linksBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "linksBed.h"

//
// Constructor
//
BedLinks::BedLinks(string &bedFile, string &base, string &org, string &db) {
    _bedFile = bedFile;
    _bed = new BedFile(bedFile);

    _base = base;
    _org = org;
    _db = db;

    CreateLinks();
}

//
// Destructor
//
BedLinks::~BedLinks(void) {
}


void BedLinks::WriteURL(BED &bed, string &base) {

    string position = bed.chrom;
    std::stringstream posStream;
    posStream << ":" << bed.start << "-" << bed.end;
    position.append(posStream.str());

    cout << "<tr>" << endl;
        cout << "\t<td>" << endl;
            cout << "\t\t<a href=" << base << position << ">";
            cout << bed.chrom << ":" << bed.start+1 << "-" << bed.end;
            cout << "</a>" << endl;
        cout << "\t</td>" << endl;

        if (_bed->bedType == 4) {
            cout << "\t<td>" << endl;
            cout << bed.name << endl;
            cout << "\t</td>" << endl;
        }
        else if (_bed->bedType == 5) {
            cout << "\t<td>" << endl;
            cout << bed.name << endl;
            cout << "\t</td>" << endl;

            cout << "\t<td>" << endl;
            cout << bed.score << endl;
            cout << "\t</td>" << endl;
        }
        else if ((_bed->bedType == 6) || (_bed->bedType == 9) || (_bed->bedType == 12)) {
            cout << "\t<td>" << endl;
            cout << bed.name << endl;
            cout << "\t</td>" << endl;

            cout << "\t<td>" << endl;
            cout << bed.score << endl;
            cout << "\t</td>" << endl;

            cout << "\t<td>" << endl;
            cout << bed.strand << endl;
            cout << "\t</td>" << endl;
        }
        cout << "</tr>" << endl;
}


void BedLinks::CreateLinks() {


    // construct the html base.
    string org = _org;
    string db = _db;
    string base = _base;
    base.append("/cgi-bin/hgTracks?org=");
    base.append(org);
    base.append("&db=");
    base.append(db);
    base.append("&position=");

    // create the HTML header
    cout << "<html>" << endl <<"\t<body>" << endl;
    cout << "<title>" << _bedFile << "</title>" << endl;

    // start the table of entries
    cout << "<br>Firefox users: Press and hold the \"apple\" or \"alt\" key and click link to open in new tab." << endl;
    cout << "<p style=\"font-family:courier\">" << endl;
    cout << "<table border=\"0\" align=\"justify\"" << endl;
    cout << "<h3>BED Entries from: stdin </h3>" << endl;


    BED bedEntry;
    _bed->Open();
    while (_bed->GetNextBed(bedEntry)) {
        if (_bed->_status == BED_VALID) {
            WriteURL(bedEntry, base);
        }
    }
    _bed->Close();

    cout << "</table>" << endl;
    cout << "</p>" << endl;
    cout << "\t</body>" << endl <<"</html>" << endl;
}


