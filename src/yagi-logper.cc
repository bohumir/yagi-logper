/*
   yagi-logper 
   Jorgen Hald, Bohumir Jelinek
   https://sourceforge.net/projects/yagi-logper/

   finds precise gain, directivity and input impedance of Yagi or
   log-periodic antennas (planar symmetric array of cylindrical
   dipoles, many driven elements) for a range of specified
   frequencies, integrates 2D for segments on the same element, moment
   method

   This software is distributed under the GNU General Public License.
   GNU GPL text is included in the LICENSE.text file.
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#ifdef _WIN32
#else
#include <X11/Xlib.h>
#include <X11/keysym.h>
#endif

using namespace std;
namespace global_names {

  static const unsigned int NDNS1MAX = 63, /* should be dynamic later */
    NDIPMAX = NDNS1MAX,      
    NSUBMAX = NDIPMAX-1,
    IW = (NSUBMAX+1)*(NSUBMAX+2),
    POINTS = 55,
    NUMFREQMAX = 51;

  static const double PI = M_PI,
    TPI = 2*M_PI,
    PIH = M_PI/2;

  typedef double my_float;
  typedef complex<my_float> my_complex;

  unsigned int ndip, nsub, ngaux, ngaufi;
  unsigned int ndns1, nexc, lexc;

  // dimensions of antenna in the multiplications of lambda
  my_float hl[NDIPMAX],hgap[NDIPMAX],posy[NDIPMAX],posz[NDIPMAX];
  my_float ar[NDIPMAX],dxk[NDIPMAX];
  unsigned int tag[NDIPMAX];
  my_complex v[NDIPMAX];
  my_float lambda;

  unsigned int numfreq,freq_cnt,pordip;
  my_float freqd,freqh,deltafreq;
  my_float dhl[NDIPMAX],dhgap[NDIPMAX],dposy[NDIPMAX],dposz[NDIPMAX];
  my_float dar[NDIPMAX];
  my_float pthetamax[2][NUMFREQMAX],pdmax[2][NUMFREQMAX];
  my_float dmax;
  my_complex pinimped[NUMFREQMAX];

  /* could be dynamic */
  float clocks_array[2][NUMFREQMAX][361];

  bool scaling;
  bool linear_gain(false);
  bool f_amplit(false), f_inimp(false), f_power(false), f_direct(false);
  bool s_amplit(false), s_inimp(false), s_power(false), s_direct(false);
  bool w_distrib(false),                                w_direct( true);

  /* could be dynamic */
  my_complex pfc[NDNS1MAX][NDNS1MAX];

  my_complex right_first[NDNS1MAX],rightz[NDNS1MAX];
  my_float amplit[NDNS1MAX],args[NDNS1MAX];
  my_float pgauss[2][16],wgauss[2][16];
  my_float power;
  my_complex inimped;

  unsigned int getmaxx = 1024;
  unsigned int getmaxy =  768;
  unsigned int nplotx = 4;
  unsigned int nploty = 3;

#ifdef _WIN32
#else
  Display* dpy;
  GC gc;
  Window w;
#endif
}

using namespace global_names;

#ifdef _WIN32
#else
/*****************************************************************************/
bool initgr() {

  dpy = XOpenDisplay(NULL);
  if (!dpy) {
    cout << "Can't open X display." << endl;
    return false;
  }
  else {

    int blackColor = BlackPixel(dpy, DefaultScreen(dpy));
    int whiteColor = WhitePixel(dpy, DefaultScreen(dpy));

    // Create the window
    //
    w = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0, 0, 
                            getmaxx, getmaxy, 0, blackColor, whiteColor);

    // Create a "Graphics Context"
    //
    gc = XCreateGC(dpy, w, 0, NULL);

    // We want to get MapNotify events
    //
    XSelectInput(dpy, w, ExposureMask | StructureNotifyMask | ButtonPressMask | KeyPressMask | KeyReleaseMask);

    // "Map" the window (that is, make it appear on the screen)
    //
    XMapWindow(dpy, w);

    return true;
  }
}

void closegraph() {
  XCloseDisplay(dpy);
}

void outtextxy(int x, int y, string str) {
  XDrawString(dpy, w, gc, x, y+10, str.c_str(), str.length());
}

void circle(int x, int y, int r) {
  XDrawArc(dpy, w, gc, x-r, y-r, 2*r, 2*r, 0, 64*360);
}

void putpixel(int x, int y, int color) {
  XDrawPoint(dpy, w, gc, x, y);
  XFlush(dpy);
}
#endif

/*****************************************************************************/

bool getTheLine(ifstream& is, string& str) {

  bool success, comment, empty;
  do {
    success = getline(is, str);
    comment = ((str.length() >= 1) && (str[0] == '#'));
    empty = (str == "");
    //    cout << "***********" << str << endl;
  } while (success && (comment || empty));
  return success;
}

bool getNECLine(ifstream& is, string& str) {

  bool success, comment, empty;
  do {
    success = getline(is, str);
    comment = (str.length() >= 2) && ((str.substr(0,2) == "CM")||(str.substr(0,2) == "CE"));

    empty = (str == "");
    //    cout << "***********" << str << endl;
  } while (success && (comment || empty));
  return success;
}

bool str2mf(const string& s, my_float& out) {
  istringstream i(s);
  i >> out;
  return !i.fail();
}

string int2string(const int& number) {

   ostringstream oss;
   oss << number;
   return oss.str();
}

string float2string(const float& number) {

   ostringstream oss;
   oss << number;
   return oss.str();
}

/* reading input file */
void input(const int argc, char* argv[] ) {

  ifstream inFile(argv[1]);
  if (!inFile) {
    cout << "Error: Input file not found !" << endl;
    exit(1);
  }

  cout << "Reading input file \"" << argv[1] << "\" ..." << endl;

  // find if the frequency is specified first
  // if not, all dimensions are assumed to be in terms of lambda
  //
  string str;
  getTheLine(inFile, str);
  string look_for = "frequency ";
  cout << str;
  if (str.find(look_for, 0) != string::npos) {
    my_float freq;
    if (!str2mf(str.substr(look_for.length()), freq)) {
      cout << "Error: Input: could not read frequency" << endl;
      cout << " from string: \"" << str << "\"" << endl;
      exit(1);
    }
    else {
      cout << "Base frequency : " << freq << " MHz" << endl;
    }

    // lambda in the untis of mm
    lambda = 299800/freq;
    getTheLine(inFile, str);
  }
  else {
    lambda = 1;
    cout << " Assuming dimensions in Lambda (frequency not specified)" << endl;
  }

  // find "dipoles" keyword with the # of dipoles followed by dimensions
  //
  look_for = "dipoles ";
  if (str.find(look_for,0) != string::npos) {
    istringstream istrm(str.substr(look_for.length()));
      
    istrm >> ndip;
    if (istrm.fail()) {
      cout << "Error: Input: could not read # of dipoles" << endl;
      cout << " from string: \"" << str << "\"" << endl;
      exit(1);
    }
    if (ndip == 0) {
      cout << "Error: Input: Number of dipoles shoudl be > 0" << endl;
      exit(1);
    }

    if (istrm.peek() == EOF) {
      nsub =  NDNS1MAX/ndip - 1;
      ngaux = 4;
      ngaufi = 4;
    }
    else {
      istrm >> nsub >> ngaux >> ngaufi;
      if (istrm.fail()) {
        cout << "Error: Input: ndip should be followed by nsub,ngaux,ngaufi";
        cout << endl;
        exit(1);
      }
    }

    ndns1 = ndip * (nsub + 1);
    if (ndns1 > NDNS1MAX) {
      cout << "Error input: ndns1 > NDNS1MAX !" << endl;
      cout << "           " << ndns1 << ">" << NDNS1MAX << endl;
      exit(1);
    }

  }
  else {
    cout << "Error: Input: unexpected input: " << endl;
    cout << " \"" << str << "\"" << endl;
    cout << " expected: dipoles (or frequency)" << endl;
    exit(1);
  }

  // read dipole positions and sizes
  //
  for (unsigned int i = 0; i < ndip; i++) {
    getTheLine(inFile, str);
    istringstream istrm2(str);

    my_float length_mm, posy_mm, posz_mm, diam_mm, gap_mm;
    
    istrm2 >> length_mm >> posy_mm >> posz_mm >> diam_mm >> gap_mm;
    istrm2 >> v[i].real() >> v[i].imag();
    if (istrm2.fail()) {
      cout << "Error: Input: could not read dipole # " << i+1 << endl;
      exit(1);
    };

    hl[i] = length_mm/2/lambda;
    posy[i] = posy_mm/lambda;
    posz[i] = posz_mm/lambda;
    ar[i] = diam_mm/2/lambda;
    hgap[i] = gap_mm/2/lambda;
  }

  // print dipoles
  //
  cout << "ndip: " << ndip << "  nsub: " << nsub <<
    "  ngaux: " << ngaux << "  ngaufi: " << ngaufi << endl;
  cout << " num  hl     posy    posz      ar        hgap      v.re   v.im"
       << endl;
  for (unsigned int i = 0; i < ndip; i++) {
    cout << fixed;
    cout << setw(3) << i + 1;
    cout << fixed;
    cout.precision(4);
    cout << " " << setw(3) << hl[i];
    cout.precision(4);
    cout << " " << setw(7) << posy[i] << " " <<  setw(7) << posz[i];
    cout << scientific;
    cout << " " << setw(7) << ar[i] << " " << hgap[i];
    cout << fixed;
    cout.precision(2);
    cout << " " << setw(6) << v[i].real();
    cout << " " << setw(6) << v[i].imag() << endl;
  }
  cout << scientific;

  // read scaling
  //
  if (!inFile.eof()) {
    getTheLine(inFile, str);
    
    string look_for = "scaling ";
    if (str.find(look_for, 0) != string::npos) {
      istringstream istrm(str.substr(look_for.length()));

      istrm >> numfreq >> freqd >> freqh;
      if (istrm.fail()) {
        cout << "Error: Input: could not read scaling numfreq,freqd,freqh from"
             << endl;
        cout << " from this string: \"" << str << "\"" << endl;
        exit(1);
      }
      cout << "Scaling numfreq: " << numfreq;
      if (numfreq >  NUMFREQMAX) {
        cout << endl;
        cout << "Error: Input: number of scaling frequencies is larger than max (";
        cout << NUMFREQMAX << ")" << endl;
        exit(1);
      }
      cout << " low: " << freqd << " hi: " << freqh << endl;
      scaling = true;
      
    }
    else {
      scaling = false;
    }
  }

  while (getTheLine(inFile, str)) {
    if (str.find("linear gain", 0) != string::npos) {
      linear_gain = true;
    }
    else if (str.find("file", 0) != string::npos) {
      if (str.find("amplit") != string::npos) {
        f_amplit = true;
      }
      if (str.find("inimp") != string::npos) {
        f_inimp = true;
      }
      if (str.find("power") != string::npos) {
        f_power = true;
      }
      if (str.find("direct") != string::npos) {
        f_direct = true;
      }
    }
    else if (str.find("screen", 0) != string::npos) {
      if (str.find("amplit") != string::npos) {
        s_amplit = true;
      }
      if (str.find("inimp") != string::npos) {
        s_inimp = true;
      }
      if (str.find("power") != string::npos) {
        s_power = true;
      }
      if (str.find("direct") != string::npos) {
        s_direct = true;
      }
    }
    else if (str.find("window", 0) != string::npos) {
      if (str.find("distrib") != string::npos) {
        w_distrib = true;
      }
      if (str.find("direct") != string::npos) {
        w_direct = true;
      }
    }
    else {
      cout << "Error: Input: unexpected input: " << endl;
      cout << " \"" << str << endl;
      cout << " expected: linear gain, file, screen or window" << endl;
      exit(1);
    }
  }

  inFile.close();

} /* input */


/* reading NEC input file */
void inputNEC(const int argc, char* argv[] ) {

  ifstream inFile(argv[1]);
  if (!inFile) {
    cout << "Error: Input: input file not found !" << endl;
    exit(1);
  }

  cout << "Reading input file \"" << argv[1] << "\" ..." << endl;

  // //find if the frequency is specified first
  // for now, all dimensions are assumed to be in terms of lambda
  //
  lambda = 1;
  cout << " Assuming dimensions in Lambda (frequency not specified)" << endl;

  string str;

  ngaux = 4;
  ngaufi = 4;

  /*
  istrm >> nsub >> ngaux >> ngaufi;
  if (istrm.fail()) {
    cout << "Error: Input: ndip should be followed by nsub,ngaux,ngaufi";
    cout << endl;
    exit(1);
  }
  */

  ndip = 0;
  bool done = false;
  my_float centy = 0;

  // read dipole positions and sizes, expect "GW" keyword
  //
  do {
    
    getNECLine(inFile, str);    

    string look_for;
    look_for = "GW ";
    if (str.find(look_for,0) != string::npos) {
      istringstream istrm(str.substr(look_for.length()));
      
      unsigned int tagc;
      istrm >> tagc;
      if (istrm.fail()) {
        cout << "Error: Input: could not read tag" << endl;
        cout << "  input line: "<< endl;
        cout << str << endl;
        exit(1);
      }
      tag[ndip] = tagc;

      unsigned int nsubc;
      istrm >> nsubc;
      if (istrm.fail()) {
        cout << "Error: Input: could not read # of segments." << endl;
        cout << "  input line: "<< endl;
        cout << str << endl;
        exit(1);
      }

      if (nsubc < 1) {
        cout << "Error: Input: nsub < 1: " << nsubc << endl;
        cout << " input line: "<< endl;
        cout << str << endl;
        exit(1);
      }

      if ((nsubc % 2) != 1) {
        cout << "Error: Input: only odd # of segments allowed." << endl
             << "  found: " << nsubc << endl;
        cout << " input line: "<< endl;
        cout << str << endl;
        exit(1);
      }

      if (ndip == 0) {
        nsub = (nsubc - 1) / 2;
      } else {
        if (((nsubc - 1)/2) != nsub) {
          cout << "Error: Input: All elements must have same # of segments."
               << endl;
          cout << "  found " << nsubc
               << "  segments, should be " << nsub
               << endl;
          cout << "  input line: "<< endl;
          cout << str << endl;
          exit(1);
        }
      }

      my_float posx_mm, posy_mm, posz_mm;
      my_float mosx_mm, mosy_mm, mosz_mm, dia_mm;
    
      istrm >> mosx_mm >> mosy_mm >> mosz_mm
            >> posx_mm >> posy_mm >> posz_mm
            >> dia_mm;
      if (istrm.fail()) {
        cout << "Error: Input: could not read pos/radius." << endl;
        cout << "  input line: "<< endl;
        cout << str << endl;
        exit(1);
      }

      if (ndip == 0) {
        centy = (posy_mm + mosy_mm)/2;
      } else {
        my_float centyc = (posy_mm + mosy_mm)/2;

        if (centyc != centy) {
          cout << "Error: Input: All elements must be symmetric wrt. xz plane!"
               << endl;
          cout << "   Element # " << ndip+1
               << " has y center " << centyc
               << ", should be " << centy
               << endl;
          exit(1);
        }

        if (posx_mm != mosx_mm) {
          cout << "Erorr: All elements must be parallel to xz plane!"
               << endl;
          cout << "   Element # " << ndip+1
               << " starts at x=" << mosx_mm << ", ends at x=" << posx_mm
               << "," << endl;
          cout << "but should start and end at the same x coord."
               << endl;
          exit(1);
        }

        if (posz_mm != mosz_mm) {
          cout << "Erorr: All elements must be parallel to xz plane!";
          cout << "   Element number " << ndip+1
               << "starts at z " << mosz_mm << " , ends at " << posz_mm
               << "but should start and end at the same z coord."
               << endl;
          exit(1);
        }
      }

      if (istrm.fail()) {
        cout << "Error: Input: could not read positions # " << ndip+1 << endl;
        cout << str << endl;
        exit(1);
      }
      //      posx[ndip] = posx_mm/lambda; must be 0
      // +/- y and z must be the same
      
      //      hl[ndip] = length_mm/2/lambda;
      hl[ndip] = (posy_mm-centy)/lambda;
      posy[ndip] = posz_mm/lambda;
      posz[ndip] = posx_mm/lambda;
      //      ar[ndip] = dia_mm/2/lambda;
      ar[ndip] = dia_mm/lambda; // nec2 expects radius
      hgap[ndip] = 0; ////gap_mm/2/lambda;
      ndip = ndip + 1;
  
    } else {
      done = true;
    }

  } while (!done);

  if (ndip == 0) {
    cout << "Error: Input: did not read any dipole." << endl;
    exit(1);
  }

  ndns1 = ndip * (nsub + 1);
  if (ndns1 > NDNS1MAX) {
    cout << "Error input: ndns1 > NDNS1MAX !" << endl;
    cout << "           " << ndns1 << ">" << NDNS1MAX << endl;
    exit(1);
  }

  string look_for;
  look_for = "GE ";
  if (str.find(look_for,0) == string::npos) {
    cout << "Error: Input: Expected \"" << look_for << "\"" << endl;
    cout << "  input line:"<< endl;
    cout << str << endl;
    exit(1);
  }

  istringstream istrm(str.substr(look_for.length()));
  unsigned int gpflag;
  istrm >> gpflag;
  if (istrm.fail()) {
    cout << "Error: Input: Could not read gpflag (unsigned int).";
    cout << "  input line:"<< endl;
    cout << " \"" << str << "\"" << endl;
    exit(1);
  }
  if (gpflag != 0) {
    cout << "Error: Input: ground plane not implemented.";
    cout << "  input line:"<< endl;
    cout << " \"" << str << "\"" << endl;
    cout << "  integer following GE should be 0" << endl;
    exit(1);
  }

  // read dipole excitation and sizes, expect "EX" keyword
  //
  done = false;
  bool source_assigned = false;
  do {

    getNECLine(inFile, str);

    look_for = "EX ";
    if (str.find(look_for,0) != string::npos) {

      istrm.str(str.substr(look_for.length()));

      unsigned int exc_type;
      istrm >> exc_type;
      if (istrm.fail()) {
        cout << "Error: Input: could not read excitation type (unsigned int).";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }
      if (exc_type != 0) {
        cout << "Error: Input: only voltage source excitation implemented.";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        cout << "  integer following EX should be 0" << endl;
      }

      unsigned int tag_number;
      istrm >> tag_number;
      if (istrm.fail()) {
        cout << "Error: Input: could not read tag number of the EX wire";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }
      if (tag_number == 0) {
        cout << "Error: No absolute segment numbering implemented.";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        cout << "  tag number should not be 0" << endl;
        exit(1);
      }

      unsigned int segment_number;
      istrm >> segment_number;
      if (istrm.fail()) {
        cout << "Error: Input: could not read segment number of the EX";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }
      if (segment_number != nsub + 1) {
        cout << "Error: Input: voltage source must be at center segment."
             << endl;
        cout << "  found: " << segment_number << endl;
        cout << "  should be: " << nsub+1 << endl;
        cout << "  input line:" << endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }
      unsigned int i4;
      istrm >> i4;
      if (istrm.fail()) {
        cout << "Error: Input: expected integer I4 field in EX";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }
      if (i4 != 0) {
        cout << "Error: Input: I4 field in EX should be 0";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }

      my_float f1, f2, f3, f4, f5, f6;
      istrm >> f1 >> f2 >> f3 >> f4 >> f5 >> f6;      
      if (istrm.fail()) {
        cout << "Error: Input: failed to read f1-f6 fields in EX";
        cout << "  input line:"<< endl;
        cout << " \"" << str << "\"" << endl;
        exit(1);
      }

      // loop over all segments looking for this tag number, set gap
      // width and voltage
      for(unsigned int cnt=0; cnt<ndip; ++cnt) {
        if (tag[cnt] == tag_number) {
          v[cnt].real() = f1;
          v[cnt].imag() = f2;
          hgap[cnt] = f6/2;
          source_assigned = true;
        }
      }
    }
    else {
      done = true;
    }
  } while (!done);

  if (source_assigned == false) {
    cout << "Error: no voltage source specified." << endl;
    exit(1);
  } else {
    cout << "Note: \"gap\" width for voltage source is the last field of EX"
         << endl
         << "  i.e. I1, I2, I3, I4, F1, F2 of the EX line are relevant"
         << endl
         << "  moreover, F6 means the width of the excitation \"gap\"."
         << endl;
  }

  // print dipoles
  //
  cout << "ndip: " << ndip << "  nsub: " << nsub <<
    "  ngaux: " << ngaux << "  ngaufi: " << ngaufi << endl;
  cout << " num  hl     posy    posz      ar        hgap      v.re   v.im"
       << endl;
  for (unsigned int i = 0; i < ndip; i++) {
    cout << fixed;
    cout << setw(3) << i + 1;
    cout << fixed;
    cout.precision(4);
    cout << " " << setw(3) << hl[i];
    cout.precision(4);
    cout << " " << setw(7) << posy[i] << " " <<  setw(7) << posz[i];
    cout << scientific;
    cout << " " << setw(7) << ar[i] << " " << hgap[i];
    cout << fixed;
    cout.precision(2);
    cout << " " << setw(6) << v[i].real();
    cout << " " << setw(6) << v[i].imag() << endl;
  }
  cout << scientific;

  // read scaling
  //
  // if (!inFile.eof()) {
  //   getTheLine(inFile, str);
    
  //   string look_for = "scaling ";
  //   if (str.find(look_for, 0) != string::npos) {
  //     istringstream istrm(str.substr(look_for.length()));

  //     istrm >> numfreq >> freqd >> freqh;
  //     if (istrm.fail()) {
  //       cout << "Error: Input: could not read scaling numfreq,freqd,freqh from"
  //            << endl;
  //       cout << " from this string: \"" << str << "\"" << endl;
  //       exit(1);
  //     }
  //     cout << "Scaling numfreq: " << numfreq;
  //     if (numfreq >  NUMFREQMAX) {
  //       cout << endl;
  //       cout << "Error: number of scaling frequencies is larger than max (";
  //       cout << NUMFREQMAX << ")" << endl;
  //       exit(1);
  //     }
  //     cout << " low: " << freqd << " hi: " << freqh << endl;
  //     scaling = true;
      
  //   }
  //   else {
  //     scaling = false;
  //   }
  // }

  // while (getTheLine(inFile, str)) {
  //   if (str.find("linear gain", 0) != string::npos) {
  //     linear_gain = true;
  //   }
  //   else if (str.find("file", 0) != string::npos) {
  //     if (str.find("amplit") != string::npos) {
  //       f_amplit = true;
  //     }
  //     if (str.find("inimp") != string::npos) {
  //       f_inimp = true;
  //     }
  //     if (str.find("power") != string::npos) {
  //       f_power = true;
  //     }
  //     if (str.find("direct") != string::npos) {
  //       f_direct = true;
  //     }
  //   }
  //   else if (str.find("screen", 0) != string::npos) {
  //     if (str.find("amplit") != string::npos) {
  //       s_amplit = true;
  //     }
  //     if (str.find("inimp") != string::npos) {
  //       s_inimp = true;
  //     }
  //     if (str.find("power") != string::npos) {
  //       s_power = true;
  //     }
  //     if (str.find("direct") != string::npos) {
  //       s_direct = true;
  //     }
  //   }
  //   else if (str.find("window", 0) != string::npos) {
  //     if (str.find("distrib") != string::npos) {
  //       w_distrib = true;
  //     }
  //     if (str.find("direct") != string::npos) {
  //       w_direct = true;
  //     }
  //   }
  //   else {
  //     cout << "Error: Input: unexpected input: " << endl;
  //     cout << " \"" << str << endl;
  //     cout << " expected: linear gain, file, screen or window" << endl;
  //     exit(1);
  //   }
  // }

  inFile.close();

} /* input */


/*****************************************************************************/
void dqgini(unsigned int * ngau, unsigned int& ngau1, unsigned int& ngau2) {
  unsigned int ngauss[2];
  
  for (unsigned int j = 0; j < 2; ++j) {
    unsigned int k = ngau[j];
    if (k == 1) {
      ngauss[j] = 1;
      pgauss[j][0] = 0.577350269189626;
      wgauss[j][0] = 1;
    }
    else if (k == 2) {
      ngauss[j] = 2;
      pgauss[j][0] = 0.339981043584856;
      pgauss[j][1] = 0.861136311594053;
      wgauss[j][0] = 0.652145154862546;
      wgauss[j][1] = 0.347854845137454;
    }
    else if (k == 3) {
      ngauss[j] = 4;
      pgauss[j][0] = 0.183434642495650;
      pgauss[j][1] = 0.525532409916329;
      pgauss[j][2] = 0.796666477413627;
      pgauss[j][3] = 0.960289856497536;
      wgauss[j][0] = 0.362683783378362;
      wgauss[j][1] = 0.313706645877887;
      wgauss[j][2] = 0.222381034453374;
      wgauss[j][3] = 0.101228536290376;
    }
    else if ((k == 4) || (k == 5) ) {
      ngauss[j] = 8;
      pgauss[j][0] = 0.9501250983763744E-1;
      pgauss[j][1] = 0.2816035507792589;
      pgauss[j][2] = 0.4580167776572274;
      pgauss[j][3] = 0.6178762444026437;
      pgauss[j][4] = 0.7554044083550030;
      pgauss[j][5] = 0.8656312023878317;
      pgauss[j][6] = 0.9445750230732326;
      pgauss[j][7] = 0.9894009349916499;
      wgauss[j][0] = 0.1894506104550685;
      wgauss[j][1] = 0.1826034150449236;
      wgauss[j][2] = 0.1691565193950025;
      wgauss[j][3] = 0.1495959888165767;
      wgauss[j][4] = 0.1246289712555339;
      wgauss[j][5] = 0.9515851168249278E-1;
      wgauss[j][6] = 0.6225352393864789E-1;
      wgauss[j][7] = 0.2715245941175409E-1;
    }
    k = ngauss[j];
    for (unsigned int i = 0; i < k; ++i) {
      pgauss[j][i] = pgauss[j][i]/2;
      wgauss[j][i] = wgauss[j][i]/2;
    }
  };
  for (unsigned int j = 0; j < 2; ++j) {
    unsigned int k = ngauss[j];
    for (unsigned int i = 0; i < k; ++i) {
      int ii = i+k;
      pgauss[j][ii] = 0.5+pgauss[j][i];
      pgauss[j][i] = 0.5-pgauss[j][i];
      wgauss[j][ii] = wgauss[j][i];
    }
  };
  ngau1 = ngauss[0]*2;
  ngau2 = ngauss[1]*2;
  for (unsigned int i = 0; i < ngau2; ++i) {
    pgauss[1][i] = sin(PIH*pgauss[1][i]);
  }
} /* dqgini */

/*****************************************************************************/
void copydip(unsigned int k, unsigned int ns1, unsigned int km) {
  unsigned int i1,i2,im,jm;
  unsigned int id;

  i1 = k*ns1;
  i2 = i1+nsub;
  id = km-k;
  id = id*ns1;

  for (unsigned int i = i1; i <= i2; ++i) {
    im = i+id;
    for (unsigned int j = i1; j <= i2; ++j) {
      jm = j+id;
      pfc[i][j] = pfc[im][jm];
    }
  }

} /* copydip */

namespace plasyd_names{
  unsigned int ngau[2];
  unsigned int ns1,ns2,ns22,i,k,ngau1,ngau2,kk,kk1,ll;
  unsigned int km;
  my_complex res;
}

namespace samdip_names {
  my_complex pw[IW];
  my_float fik,fi,sfik,alfak,beta2,bs,bls,s,t,tt,q;
  my_complex p,c,cc,x,y,z;
  unsigned int i,n,n1,ns21,i0,ii,jj,j;
}

void f00(int i) {
  using namespace plasyd_names;
  using namespace samdip_names;
  my_float x,xx,a,b,r;

  res.real() = 0;
  res.imag() = 0;
  x = pgauss[0][i];
  xx = x*x;
  a = sin(fik*(1-x))/sfik;
  b = 1-q*x;
  for (unsigned int j = 0; j < ngau2; ++j) {
    s = alfak*pgauss[1][j];
    r = sqrt(xx+s*s);
    t = fik*r;
    if (t>1E-6) {
      s = 1/r;
      res.real() = res.real()+(a*s*cos(t)-b*s)*wgauss[1][j];
      res.imag() = res.imag()+(-a*s*sin(t)+b*fik)*wgauss[1][j];
    }
  }
} /* f00 */

void f01a(int i) {
  using namespace plasyd_names;
  using namespace samdip_names;
  my_float x,xx,a,r,t;

  res.real() = 0;
  res.imag() = 0;
  x = pgauss[0][i];
  xx = (x+1)*(x+1);
  a = sin(fik*(1-x))/sfik;
  for (unsigned int j = 0; j < ngau2; ++j) {
    s = alfak*pgauss[1][j];
    r = sqrt(xx+s*s);
    t = fik*r;
    res.real() = res.real()+cos(t)*wgauss[1][j]/r;
    res.imag() = res.imag()-sin(t)*wgauss[1][j]/r;
  };
  res.real() = res.real()*a;
  res.imag() = res.imag()*a;
} /* f01a */

void f01b(int i) {
  using namespace plasyd_names;
  using namespace samdip_names;
  my_float x,xx,a,b,r,t;

  res.real() = 0;
  res.imag() = 0;
  x = pgauss[0][i];
  xx = x*x;
  b = fik*x;
  a = sin(b)/sfik;
  b = b/sfik;
  for (unsigned int j = 0; j < ngau2; ++j) {
    s = alfak*pgauss[1][j];
    r = sqrt(xx+s*s);
    t = fik*r;
    if (t>1E-6) {
      res.real() = res.real()+(a*cos(t)-b)*wgauss[1][j]/r;
      res.imag() = res.imag()+(-a*sin(t)+b*t)*wgauss[1][j]/r;
    }
  }
} /* f01b */

void f0n(int i) {
  using namespace plasyd_names;
  using namespace samdip_names;

  my_float x,xxm,xxp,a,s,s2,rp,rm,tp,tm;

  res.real() = 0;
  res.imag() = 0;
  x = pgauss[0][i];
  xxm = (x-q)*(x-q);
  xxp = (x+q)*(x+q);
  a = sin(fik*(1-x))/sfik;
  for (unsigned int j = 0; j < ngau2; ++j) {
    s = alfak*pgauss[1][j];
    s2 = s*s;
    rp = sqrt(xxp+s2);
    rm = sqrt(xxm+s2);
    tp = fik*rp;
    tm = fik*rm;
    res.real() = res.real()+( cos(tp)/rp +cos(tm)/rm )*wgauss[1][j];
    res.imag() = res.imag()+(-sin(tp)/rp -sin(tm)/rm )*wgauss[1][j];
  }
  res.real() = res.real()*a;
  res.imag() = res.imag()*a;
} /* f0n */

/*****************************************************************************/
void samdip(unsigned int k) {

  using namespace plasyd_names;
  using namespace samdip_names;

  //  new(pw);
  fik = dxk[k];
  fi = fik;
  sfik = sin(fik);
  alfak = 8*PIH*ar[k]/fik;
  beta2 = alfak*alfak/2;

  bs = 0;
  bls = 0;
  for (unsigned int i = 0; i < ngau2; ++i) {
    s = alfak*pgauss[1][i];
    t = sqrt(1+s*s);
    tt = log(1+t);
    bs = bs+t*wgauss[1][i];
    bls = bls+tt*wgauss[1][i];
  };
  bs = bs-alfak/PIH; /* exactly K1,K2 */

  /* self coupling */
  q = fik*cos(fik)/sin(fik);
  s = log(2/alfak) + bls - q*bs;
  t = fik*(q/2 - 1);
  p.real() = s;
  p.imag() = t;
  c = p;

  /* dqg(@f00,c); */
  for (unsigned int i = 0; i < ngau1; ++i) {
    f00(i);
    c.real() = c.real()+wgauss[0][i]*res.real();
    c.imag() = c.imag()+wgauss[0][i]*res.imag();
  };
  pw[0].real() = 2*c.real();
  pw[0].imag() = 2*c.imag();
  /* coupling between neigbour exp.f. */
  p.real() = bs*fik/sfik;
  p.imag() = -fik/2*fik/sfik;

  /* dqg(@f01a,c); */
  my_complex c;
  c.real() = 0;
  c.imag() = 0;
  for (unsigned int i = 0; i < ngau1; ++i) {
      f01a(i);
      c.real() = c.real()+wgauss[0][i]*res.real();
      c.imag() = c.imag()+wgauss[0][i]*res.imag();
  }

  /* dqg(@f01b,cc); */
  cc = p;
  for (unsigned int i = 0; i < ngau1; ++i) {
      f01b(i);
      cc.real() = cc.real()+wgauss[0][i]*res.real();
      cc.imag() = cc.imag()+wgauss[0][i]*res.imag();
  }

  pw[1].real() = c.real()+cc.real();
  pw[1].imag() = c.imag()+cc.imag();

  /* remaining coupling integrals */
  for (unsigned int n1 = 2; n1 < ns22; ++n1) {
    n = n1-1;
    q = n+1;
    /* dqg(@f0n,c); */
    my_complex c;
    c.real() = 0;
    c.imag() = 0;
    for (unsigned int i = 0; i < ngau1; ++i) {
      f0n(i);
      c.real() = c.real()+wgauss[0][i]*res.real();
      c.imag() = c.imag()+wgauss[0][i]*res.imag();
    };
    pw[n1] = c;
  }

  /* mutual impedances from coupling integrals */
  x = pw[1];
  y = pw[0];
  t = 2*cos(fik);
  ns21 = ns22-1;
  for (unsigned int i = 0; i < ns21; ++i) {
    z = pw[i+1];
    pw[i].real() = x.real()+z.real()-t*y.real();
    pw[i].imag() = x.imag()+z.imag()-t*y.imag();
    x = y;
    y = z;
  };

  /* utilize antenna symmetry while forming impedance submatrix */
  i0 = k*ns1;
  for (unsigned int ii = 0; ii < ns1; ++ii) {
    int i = i0+ii;
    pfc[i][i].real() = pw[0].real()+pw[ii+ii].real();
    pfc[i][i].imag() = pw[0].imag()+pw[ii+ii].imag();
  }

  for (unsigned int ii = 0; ii < nsub; ++ii) {
    int i = i0+ii+1;
    for (jj = 0; jj <= ii; ++jj) {
      j = i0+jj;
      pfc[i][j].real() = pw[ii-jj+1].real()+pw[ii+jj+1].real();
      pfc[i][j].imag() = pw[ii-jj+1].imag()+pw[ii+jj+1].imag();
      pfc[j][i].real() = pfc[i][j].real();
      pfc[j][i].imag() = pfc[i][j].imag();
    }
  }
  // dispose(pw);
} /* samdip */

namespace difdip_names {
  unsigned int k,l,i,j,i0,j0,jj,m,ii,iip,iim,is,js;
  my_float ss,fi,sfi,beta2,del,q,x0,x01,dkl,t;
  my_complex c;
  my_complex pw[NSUBMAX+2][NSUBMAX+1];
}

void edif4(int m) {
  using namespace plasyd_names;
  using namespace difdip_names;
  my_float x,x1,x2,x3,x4,r1,r2,r3,r4;

  x = pgauss[0][m];
  x1 = x+x0;
  x2 = x-x0;
  x3 = x+x01;
  x4 = x-x01;
  r1 = sqrt(x1*x1+beta2);
  r2 = sqrt(x2*x2+beta2);
  r3 = sqrt(x3*x3+beta2);
  r4 = sqrt(x4*x4+beta2);
  x1 = fi*r1;   /* into numerator */
  x2 = fi*r2;
  x3 = fi*r3;
  x4 = fi*r4;
  res.real() = sin(fi*(1-x))*(cos(x1)/r1+cos(x2)/r2
                              +cos(x3)/r3+cos(x4)/r4);
  res.imag() = -sin(fi*(1-x))*(sin(x1)/r1+sin(x2)/r2
                               +sin(x3)/r3+sin(x4)/r4);
} /* edif4 */

void difdip(unsigned int kk, unsigned int ll) {
  using namespace plasyd_names;
  using namespace difdip_names;

  //  label zac,kon;
  //  new(pw);
  unsigned int k = kk;
  unsigned int l = ll;
  if (k == l) {
    cout << "Error: k == l in difdip" << endl;
    exit(1);
  }
  else {
    ss = (posy[k]-posy[l])*(posy[k]-posy[l])+
      (posz[k]-posz[l])*(posz[k]-posz[l])+
      ar[k]*ar[k]+ar[l]*ar[l];
    ss = ss*TPI*TPI;
  zac: 
    fi = dxk[l]; /* 2*pi*delx[l] */
    sfi = sin(fi);
    beta2 = ss/fi/fi;
    dkl = dxk[k]/fi;
    /* coupling integral by numerical integration */
    for (unsigned int i = 0; i < ns2; ++i) {        /* m = relatively to*/
      del = i*dkl;
      for (unsigned int j = 0; j < ns1; ++j) {      /* i = from pair*/
        q = j;
        x0 = q-del;
        x01 = q+del;
        /* dqg(edif4,c); */
        c.real() = 0;
        c.imag() = 0;
        for (unsigned int m = 0; m < ngau1; ++m) {
          edif4(m);
          c.real() = c.real()+wgauss[0][m]*res.real();
          c.imag() = c.imag()+wgauss[0][m]*res.imag();
        }
        pw[i][j].real() = c.real()/sfi;
        pw[i][j].imag() = c.imag()/sfi;
      }
    }

    /* use integrals to form impedance matrix elements */
    unsigned int i0 = k*ns1;
    unsigned int j0 = l*ns1;
    unsigned int i = i0;
    t = cos(dxk[k]);
    for (unsigned int jj = 0; jj < ns1; ++jj) {
      j = j0+jj;
      pfc[i][j].real() = 2*(pw[1][jj].real()-t*pw[0][jj].real());
      pfc[i][j].imag() = 2*(pw[1][jj].imag()-t*pw[0][jj].imag());
    };
    t = 2*t;
    for (unsigned int ii = 1; ii < ns1; ++ii) {
      i = i0+ii;
      iip = ii+1;
      iim = ii-1;
      for (unsigned int jj = 0; jj < ns1; ++jj) {
        j = j0+jj;
        pfc[i][j].real() = pw[iip][jj].real()+pw[iim][jj].real()-t*pw[ii][jj].real();
        pfc[i][j].imag() = pw[iip][jj].imag()+pw[iim][jj].imag()-t*pw[ii][jj].imag();
      }
    }

    /* interchange k and l */
    if (k == ll) {
      /* nothing */
    }
    else {
      if (dkl == 1) {
        goto kon;
      }
      else {
        k = ll;
        l = kk;
        goto zac;
      };
    kon: /* sym */
      for (unsigned int ii = 0; ii < ns1; ++ii) {
        i = i0+ii;
        is = j0+ii;
        for (unsigned int jj = 0; jj < ns1; ++jj) {
          j = j0+jj;
          js = i0+jj;
          pfc[is][js] = pfc[i][j];
        }
      }
    }
  }
  // dispose(pw);
} /* difdip */

/*****************************************************************************/
void plasyd() { /* calculates impedance matrix */

  using namespace plasyd_names;
  
  ////  cout << "Calculating impedance matrix." << endl;
  ns1 = nsub+1;
  ns2 = ns1+1;
  ns22 = ns1+ns1;
  for (unsigned int i = 0; i < ndip; ++i) {
    dxk[i] = hl[i]*TPI/ns1;
  }
  ngau[0] = ngaux+1;
  ngau[1] = ngaufi;
  dqgini(ngau, ngau1, ngau2);

  /* same dipole */
  for (unsigned int k = 0; k < ndip; ++k) {
    unsigned int i = 0;
    while (i<k) {
      if ((dxk[k] == dxk[i]) && (ar[i] == ar[k])) {
        km = i;
        copydip(k,ns1,km);
        i = k;
      };
      ++i;
    }
    if (i == k) {
      samdip(k);
    }
  }

  /* different dipoles */
  if (ndip != 1) {
    ngau[0] = ngaux;
    ngau[1] = 1;
    dqgini(ngau,ngau1,ngau2);
    for (unsigned int kk = 1; kk < ndip; ++kk) {
      for (unsigned int ll = 0; ll < kk; ++ll) {
        difdip(kk,ll);
      }
    }
  }

} /* plasyd */

/*****************************************************************************/
namespace cgapex_names {
  unsigned int i,k,m,ns1,kr;
  my_float fik,gk,g0k,x1mg,x2mg,x3mg;
  my_complex c;
}

void shift() {
  using namespace cgapex_names;
  m = m+1;
  x1mg = x2mg;
  x2mg = x3mg;
  x3mg = x3mg+fik;
}

void cgapex() {
  using namespace cgapex_names;
  for (unsigned int i = 0; i < ndns1; ++i) {
    right_first[i].real() = 0;
    right_first[i].imag() = 0;
  };
  unsigned int k = 0;
  ns1 = nsub+1;
  nexc = 0;
  /* dip. loop starts here */
  do {
    if ((v[k].real() != 0) || (v[k].imag() != 0)) {
      nexc = nexc+1;
      lexc = k;
      fik = dxk[k];
      gk = TPI*hgap[k];
      c.real() = v[k].imag()/gk/30;
      c.imag() = -v[k].real()/gk/30;
      g0k = gk;
      if (fik<gk) {
        g0k = fik;
      }
      kr = k*ns1;
      right_first[kr].real() = 2*c.real()*sin(g0k/2)*sin(fik-g0k/2);
      right_first[kr].imag() = 2*c.imag()*sin(g0k/2)*sin(fik-g0k/2);
      m = 0;
      x1mg = -gk;
      x2mg = x1mg+fik;
      x3mg = x2mg+fik;
      while (m < nsub) {
        kr = kr+1;
        if (x3mg <= 0 ) {
          right_first[kr].real() = right_first[kr-1].real();
          right_first[kr].imag() = right_first[kr-1].imag();
          shift();
        }
        else {
          if (x2mg <= 0) {
            right_first[kr].real() = right_first[kr-1].real()/2+
              c.real()*sin(-x2mg/2)*sin((x3mg+fik)/2);
            right_first[kr].imag() = right_first[kr-1].imag()/2+
              c.imag()*sin(-x2mg/2)*sin((x3mg+fik)/2);
            shift();
          }
          else {
            if (x1mg < 0) {
              right_first[kr].real() = c.real()*sin(-x1mg/2)*
                sin(-x1mg/2);
              right_first[kr].imag() = c.imag()*sin(-x1mg/2)*
                sin(-x1mg/2);
              shift();
            }
            else {
              m = nsub+1;
            }
          }
        }
      }
    }
    k = k + 1;
  } while (k < ndip);
  if (nexc == 0) {
    cout << "Error: No excitation on dipoles !" << endl;
    exit(1);
  }
} /* cgapex */
  
/*****************************************************************************/
void cgauss() {
  my_float h;
  my_complex c;

  ////  cout << "Solving the system of equations." << endl;

  for (unsigned int jz = 0; jz < ndns1; ++jz) {
    h = pfc[jz][jz].real()*pfc[jz][jz].real()+pfc[jz][jz].imag()*pfc[jz][jz].imag();
    c.real() = (right_first[jz].real()*pfc[jz][jz].real()
                +right_first[jz].imag()*pfc[jz][jz].imag())/h;
    c.imag() = (right_first[jz].imag()*pfc[jz][jz].real()-
                right_first[jz].real()*pfc[jz][jz].imag())/h;
    right_first[jz] = c;
    for (unsigned int jh = ndns1; jh > jz; --jh) {
      c.real() = (pfc[jz][jh-1].real()*pfc[jz][jz].real()
                  +pfc[jz][jh-1].imag()*pfc[jz][jz].imag())/h;
      c.imag() = (pfc[jz][jh-1].imag()*pfc[jz][jz].real()-
                  pfc[jz][jh-1].real()*pfc[jz][jz].imag())/h;
      pfc[jz][jh-1] = c;
    }
    for (unsigned int iz = jz+1; iz < ndns1; ++iz) {
      /* right_first[iz] = right_first[iz]-right_first[jz]*pfc[iz][jz]; */
      right_first[iz].real() = right_first[iz].real()-right_first[jz].real()*pfc[iz][jz].real()
        +right_first[jz].imag()*pfc[iz][jz].imag();
      right_first[iz].imag() = right_first[iz].imag()-right_first[jz].real()*pfc[iz][jz].imag()
        -right_first[jz].imag()*pfc[iz][jz].real();
      for (unsigned int jh = ndns1-1; jh > jz; --jh) {
        pfc[iz][jh].real() = pfc[iz][jh].real()-pfc[jz][jh].real()*pfc[iz][jz].real()
          +pfc[jz][jh].imag()*pfc[iz][jz].imag();
        pfc[iz][jh].imag() = pfc[iz][jh].imag()-pfc[jz][jh].real()*pfc[iz][jz].imag()
          -pfc[jz][jh].imag()*pfc[iz][jz].real();
      };
      unsigned int jh = jz;
      c.real() = -pfc[jz][jh].real()*pfc[iz][jz].real()
        +pfc[jz][jh].imag()*pfc[iz][jz].imag();
      c.imag() = -pfc[jz][jh].real()*pfc[iz][jz].imag()
        -pfc[jz][jh].imag()*pfc[iz][jz].real();
      pfc[iz][jh].real() = pfc[iz][jh].real()+c.real();
      pfc[iz][jh].imag() = pfc[iz][jh].imag()+c.imag();
    };
  };

  /* zero out above the diagonal */
  for (unsigned int jz = ndns1-1; jz > 0; --jz) {
    c.real() = right_first[jz-1].real();
    c.imag() = right_first[jz-1].imag();
    for (unsigned int jh = ndns1-1; jh >= jz; --jh) {
      /* c = c-right_first[jh]*pfc[jz][jh]; */
      c.real() = c.real()-right_first[jh].real()*pfc[jz-1][jh].real()
        +right_first[jh].imag()*pfc[jz-1][jh].imag();
      c.imag() = c.imag()-right_first[jh].real()*pfc[jz-1][jh].imag()
        -right_first[jh].imag()*pfc[jz-1][jh].real();
    };
    right_first[jz-1] = c;
    c.real() = 0;
    c.imag() = 0;
    for (unsigned int jh = ndns1-1; jh >= jz; --jh) {
      pfc[jz-1][jh] = c;
    }
  }
} /* cgauss */

/*****************************************************************************/
void camplit() { /* calculates amplitudes and phases of currents */
  my_float h;

  ////  cout << "Calculating amplitudes." << endl;

  for (unsigned int i = 0; i < ndns1; ++i) {

    h = right_first[i].real()*right_first[i].real()+right_first[i].imag()*right_first[i].imag();
    h = sqrt(h);
    args[i] = atan(right_first[i].imag()/right_first[i].real())*180/PI;
    if (right_first[i].real() < 0) {
      args[i] = 180+args[i];
    }
    else {
      if (right_first[i].imag() < 0) {
        args[i] = 360+args[i];
      }
    }

    string cnt;
    if ((i) % (nsub+1) == 0) {
      amplit[i] = 2*h;
    }
    else {
      amplit[i] = h;
    };
    
    if (s_amplit) {
      cout << "Ampl." << i << ".: " << setprecision(15) << setw(22) << amplit[i];
      cout << "  Ph.: " << setw(22) << args[i];
      if ((i) % (nsub+1) == 0) {
        cout << " central";
      }
      cout << endl;
    }
  }
} /* camplit */

/*****************************************************************************/
/*
procedure graphdistrib;
    var i,j,ir,iz,pol: byte;
        xp,xm,ipz,xz,yz: word;
        xs,ys,dx,fik,
        sfik,ddf,ddx,xm0,xp0: single;
        st: string;
        p: array[NDNS1MAX][POINTS+1] of single;
        s: array[POINTS] of single;
        k,kp: my_type;

      if (initgr) {
        
          xs = (getmaxx+1)/2;
          dx = xs/(nsub+1);
          ys = (getmaxy+1)/2;
          for j = 1 to ndip do
              {
                str(j,st);
                setcolor(white);
                outtextxy(10,10,"Dipole number : "+st);
                line(0,round(ys),getmaxx,round(ys));
                iz = (j-1)*(nsub+1);
                k = amplit[iz+1];
                kp = args[iz+1];
                str(k:14,st);
                outtextxy(10,20,"Centr curr amplit :  "+st);
                str(round(kp*1E5)/1E5:9:5,st);
                outtextxy(10,30,"            phase : "+st);
                if ((nexc == 1) && (lexc == j)
                  ) {
                       
                         str(inimped.real():14,st);
                         outtextxy(10,40,"Input impedance (re):" +st);
                         str(inimped.imag():14,st);
                         outtextxy(10,50,"                (im):" +st)
                       };
                for i = 2 to nsub+1 do
                    {
                      ir = i+iz;
                      if (amplit[ir]>k
                        ) {
                           k = amplit[ir];
                      if (args[ir]>kp
                        ) {
                           kp = args[ir]
                    };

                str(k,st);
                outtextxy(10,40,"Curr - max amplit : "+st);
                str(kp,st);
                outtextxy(10,50,"            phase : "+st);

                fik = dxk[j];
                sfik = sin(fik);
                ddf = fik/(POINTS+1);
                ipz = iz*(POINTS+1)+1;
                p[iz+1][1] = amplit[iz+1];
                for i = 1 to POINTS do
                    {
                      s[i] = sin(fik-i*ddf)/sfik;
                      p[iz+1][i+1] = s[i]*amplit[iz+1]
                    };

                k = ys/k/2;
                kp = ys/kp/2;
                for i = 2 to nsub+1 do
                    {
                      p[iz+i][1] = amplit[iz+i];
                      for pol = 2 to POINTS+1 do
                          {
                            p[iz+i][pol] = s[pol-1]*amplit[iz+i];
                            p[iz+i-1][POINTS+3-pol] = p[iz+i-1][POINTS+3-pol]
                                                    +s[pol-1]*amplit[iz+i]
                          }
                    };
                ddx = dx/(points+1);
                for i = 1 to nsub+1 do
                   {
                      xp0 = xs+(i-1)*dx;
                      xm0 = xs-(i-1)*dx;
                      xp = round(xp0);
                      xm = round(xm0);
                      putpixel(xp,round(ys),brown);
                      putpixel(xm,round(ys),brown);
                      ir = i+iz;
                      if (i == 1) {
                             
                               right_first[ir].real() = right_first[ir].real()*2;
                               right_first[ir].imag() = right_first[ir].imag()*2
                             };
                      setcolor(blue);
                      circle(xp,round(ys-k*right_first[ir].real()),2);
                      circle(xm,round(ys-k*right_first[ir].real()),2);
                      setcolor(red);
                      circle(xp,round(ys-k*right_first[ir].imag()),2);
                      circle(xm,round(ys-k*right_first[ir].imag()),2);
                      if (i==1) {
                               right_first[ir].real() = right_first[ir].real()/2;
                               right_first[ir].imag() = right_first[ir].imag()/2
                             };
                      setcolor(white);
                      circle(xp,round(ys-k*amplit[ir]),3);
                      circle(xm,round(ys-k*amplit[ir]),3);
                      for pol = 2 to POINTS+1 do
                          {
                            xz = round(xp0+(pol-1)*ddx);
                            yz = round(ys-k*p[ir,pol]);
                            putpixel(xz,yz,magenta);
                            putpixel(round(xm0-(pol-1)*ddx),round(ys-k*p[ir,pol]),magenta)
                          };
                      setcolor(green);
                      circle(xp,round(ys-kp*args[ir]),1);
                      circle(xm,round(ys-kp*args[ir]),1)
                   };
                readkey;
                cleardevice;
              };
          CloseGraph;
        }
};*/ /*graphdistrib*/

/*****************************************************************************/
void inimp() {
  my_float pom;
  unsigned int iz;
  if (nexc == 1) {
    iz = (nsub+1)*lexc;
    pom = right_first[iz].real()*right_first[iz].real()
      +right_first[iz].imag()*right_first[iz].imag();
    inimped.real() = (v[lexc].real()*right_first[iz].real()
                      +v[lexc].imag()*right_first[iz].imag())/pom/2;
    inimped.imag() = (v[lexc].imag()*right_first[iz].real()
                      -v[lexc].real()*right_first[iz].imag())/pom/2;
    if (s_inimp) {
      cout << endl;
      cout << "Input impedance (re): " << inimped.real() << endl;
      cout << "                (im): " << inimped.imag() << endl;
    };
    if (scaling) {
      pinimped[freq_cnt] = inimped;
    }
  }
}

/*****************************************************************************/
void pabs() { /* powers */
  my_float w[NDIPMAX];
  my_float cim;
  unsigned int ns1,kr;

  power = 0;
  ns1 = nsub+1;
  for (unsigned int k = 0; k < ndip; ++k) {
    kr = k*ns1;
    cim = -right_first[kr].imag()*rightz[kr].real()
      +right_first[kr].real()*rightz[kr].imag();
    for (unsigned int m = 0; m < nsub; ++m) {
      kr = kr+1;
      cim = cim-right_first[kr].imag()*rightz[kr].real()
        +right_first[kr].real()*rightz[kr].imag();
    };
    w[k] = -30*cim/sin(dxk[k]);
    if (s_power) {
      cout << " " << k << ". dip. power: " << w[k] << endl;
    };
    power = power+w[k];
  };
  if (s_power) {
    cout << "   Total power: " << power << endl;
  };
} /* pabs */

#ifdef _WIN32
#else
/*****************************************************************************/
void polar(my_float fi, float* hodndir,
           unsigned int xs, unsigned int ys, unsigned int r, my_float dmax) {

  const int NKRU = 5;
  int dbkru[NKRU] = { 3, 10, 20, 30, 40};
  int dakru[NKRU] = { 3, 3, 6, 6, 15};
  int alfa, i, j, x, y, c;
  float alfar, rm, rca, rsa, rmm, koef, expon;
  string st;

  if (linear_gain) {
    expon = 1;
  }
  else {
    expon = 1./15*10/log(10);
  }
  koef = expon*log(10)/10;

  // setcolor(lightgreen);
  // :7:4
  st = float2string(10.*log(dmax)/log(10));
  outtextxy(xs-r,ys-r,"Max. gain: "+st+"dB");

  if (fi == 0) {
    outtextxy(xs-r,ys-r+10,"E-plane");
  }
  else {
    outtextxy(xs-r,ys-r+10,"H-plane");
  }

  if (scaling) {
    // 4:2
    if (lambda != 1) {
      st = float2string((freqd + freq_cnt*deltafreq)*299800/lambda);
      outtextxy(xs-r,ys-r+20,st+"MHz");
    }
    else {
      st = float2string(freqd + freq_cnt*deltafreq);
      outtextxy(xs-r,ys-r+20,"Sc: "+st);
    }
  };

  // setcolor(white);
  circle(xs,ys,r);
  int white = 0;
  putpixel(xs,ys,white);
  for (c = 0; c < NKRU; ++c) {
    rm = exp(-koef*dbkru[c]);
    for (i = 0; i < (45. / dakru[c]); ++i) {
      alfa = i*dakru[c];
      alfar = alfa/180.*PI;
      rca = r*cos(alfar);
      rsa = r*sin(alfar);
      if ((alfa % 30 == 0) && (c == 0)) {
        for(j = 0; j < 30; ++j) {
          rmm = exp(-koef*j);
          x = static_cast<int>(round(rmm*rca));
          y = static_cast<int>(round(rmm*rsa));
          putpixel(xs+x,ys+y,white);
          putpixel(xs+x,ys-y,white);
          putpixel(xs-x,ys+y,white);
          putpixel(xs-x,ys-y,white);
          putpixel(xs+y,ys+x,white);
          putpixel(xs+y,ys-x,white);
          putpixel(xs-y,ys+x,white);
          putpixel(xs-y,ys-x,white);
        }
      }
      x = static_cast<int>(round(rm*rca));
      y = static_cast<int>(round(rm*rsa));
      putpixel(xs+x,ys+y,white);
      putpixel(xs+x,ys-y,white);
      putpixel(xs-x,ys+y,white);
      putpixel(xs-x,ys-y,white);
      putpixel(xs+y,ys+x,white);
      putpixel(xs+y,ys-x,white);
      putpixel(xs-y,ys+x,white);
      putpixel(xs-y,ys-x,white);
    }
  }

  for (alfa = 0; alfa <= 360; ++alfa) {
    alfar = alfa/180.*PI;
    if (linear_gain) {
      rm = hodndir[alfa]/dmax;
    }
    else {
      rm = exp(expon*log(hodndir[alfa]/dmax));
    }
    //      setcolor(yellow);
    if (not scaling) {
      circle(xs+static_cast<int>(round(r*rm*cos(alfar))),
             ys+static_cast<int>(round(r*rm*sin(alfar))),
             1);
    }
    else {
      // yellow
      putpixel(xs+static_cast<int>(round(r*rm*cos(alfar))),
               ys+static_cast<int>(round(r*rm*sin(alfar))),0);
      // setcolor(white);
    }
  }
} /*polar*/
#endif

namespace gc_names {
  const int STEP = 1;
  my_float theta,thetamax,fi,delth,d;
  unsigned int np ,i;
  short int pom;
}

void direct() {
  using namespace gc_names;

  my_complex c,ck;
  my_float thr,fir,cth,sth,cfi,sfi,stsf,stcf,ap,am,s,fik;
  short int kr,ns1;

  thr = TPI*theta/360;
  fir = TPI*fi/360;
  cth = cos(thr);
  sth = sin(thr);
  cfi = cos(fir);
  sfi = sin(fir);
  stsf = sth*sfi;
  stcf = sth*cfi;
  ap = (1+stcf)/2;
  am = (1-stcf)/2;
  c.real() = 0;
  c.imag() = 0;
  ns1 = nsub+1;

  for (unsigned int k = 0; k < ndip; ++k) {
    fik = dxk[k];
    kr = k*ns1;
    s = stcf*fik;
    ck = right_first[kr];
    for (unsigned int m = 1; m <= nsub; ++m) {
      ck.real() = ck.real()+right_first[kr+m].real()*cos(m*s);
      ck.imag() = ck.imag()+right_first[kr+m].imag()*cos(m*s);
    }
    s = fik*ap;
    if (s<1E-4) {
      ck.real() = ck.real()*(1-s*s/6);
      ck.imag() = ck.imag()*(1-s*s/6);
    }
    else {
      ck.real() = ck.real()*sin(s)/s;
      ck.imag() = ck.imag()*sin(s)/s;
    };
    s = fik*am;
    if (s<1E-4) {
      ck.real() = ck.real()*(1-s*s/6);
      ck.imag() = ck.imag()*(1-s*s/6);
    }
    else {
      ck.real() = ck.real()*sin(s)/s;
      ck.imag() = ck.imag()*sin(s)/s;
    };
    ck.real() = ck.real()*fik*fik/sin(fik);
    ck.imag() = ck.imag()*fik*fik/sin(fik);
    s = TPI*(posy[k]*stsf+posz[k]*cth);
    c.real() = c.real()+ck.real()*cos(s)-ck.imag()*sin(s);
    c.imag() = c.imag()+ck.real()*sin(s)+ck.imag()*cos(s);
  };
  d = 60*(cth*cth*cfi*cfi+sfi*sfi)*(c.real()*c.real()+c.imag()*c.imag())/power;
  if (d<1E-20) {
    d = 1E-20;
  }; /* direct */

}

/*****************************************************************************/
void get_clock(char* argv[], int freq_cnt) {

  using namespace gc_names;

  fi = 0;
  np = 360 / STEP + 1;
  delth = 360 / (np-1);

  for (short int pom = 0; pom < 2; ++pom) {
    dmax = 0;
    for (unsigned int i = 0; i < np; ++i) {

      theta = i*delth;
      direct();

      clocks_array[pom][freq_cnt][i] = d;

      if (d>dmax) {
        dmax = d;
        thetamax = theta;
      };

    }

    if (scaling) {
      if (thetamax == 360.0) {
        thetamax = 0;
      }
      pthetamax[pom][freq_cnt] = thetamax;
      pdmax[pom][freq_cnt] = dmax;
    }
    else {
      if (! scaling) {
#ifdef _WIN32
#else
        if (initgr()) {

          // Wait for the MapNotify event
          //
          for(;;) {
            XEvent e;
            XNextEvent(dpy, &e);
            //cout << "got next event, looking for Expose: " << Expose << endl;
            //cout << e.type << endl;
            if ((e.type == Expose) && e.xexpose.count == 0) {
              //cout << getmaxx / 4;
              polar(fi,clocks_array[pom][0],
                    getmaxx / 2,
                    getmaxy / 2,
                    static_cast<int>(round(getmaxy/2.)),
                    dmax);
              //      cout << "here" << endl;
              //XDrawLine(dpy, w, gc, 10, 60, 180, 20);
              //XDrawArc(dpy, w, gc, 5, 55, 10, 10, 0, 48*360);
              //XDrawString(dpy, w, gc, 2, 20, "Hello World!", 12);
              //XFlush(dpy);
            }
            else if(e.type==ButtonPress) {
              //cout << "Buttonpress found" << endl;
              break;
            }
            else if (e.type==KeyPress) {
              //cout << "keyPress found" << endl;
              //cout << e.xkey.state << endl;
              //cout << XKeycodeToKeysym(dpy,e.xkey.keycode,0) << endl;
              if (XLookupKeysym(&e.xkey, 0) != XK_Alt_L) {
                break;
              }
            }
          }
          closegraph();
        }
#endif
      }
    }
    fi = 90;
  }

  if (f_direct) {
    
    for (short int pom = 0; pom < 2; ++pom) {

      string ofn(argv[1]);
      if (linear_gain) {
        ofn += ".dirlin.";
      }
      else {
        ofn += ".dirlog.";
      }

      if (pom == 0) {
        ofn += "E";
      }
      else {
        ofn += "H";
      }

      if (scaling) {
        ofn += ".";
        ofn += int2string(freq_cnt);
      }

      ofstream outFile(ofn.c_str());

      for (unsigned int i = 0; i < np; ++i) {

        theta = i*delth;
        my_float d = clocks_array[pom][freq_cnt][i];

        outFile << fixed << setprecision(1);
        outFile << setw(5) << theta;
        if (linear_gain) {
          outFile << setprecision(6) << setw(11) << d
                  << endl;
        }
        else {
          outFile << setprecision(6) << setw(10) << " " << 10*log(d)/log(10)
                  << endl;
        }
      }
      outFile.close();
    }
  }

} /* get_clock */

#ifdef _WIN32
#else
/*****************************************************************************/
void manyplots() {
  unsigned int r, i, j, k;

  if (initgr()) {

    r = getmaxx / (2*nplotx);
    if ((getmaxy / (2*nploty)) < r) {
      r = getmaxy / nploty;
    }
  
    // E and H plane
    bool newdraw = false;
    for (k = 0; k < 2; ++k) {

      freq_cnt = 0;
      int freq_curr_start = 0;
      while (freq_cnt < numfreq) {

        while(1) {
          XEvent e;

          if (! newdraw) {
            XNextEvent(dpy, &e);
            // cout << "got next event, looking for Expose: " << Expose << endl;
            // cout << e.type << endl;
          }

          if (newdraw || ((e.type == Expose) && (e.xexpose.count == 0))) {

            if (!newdraw) { // i.e. only expose
              freq_cnt = freq_curr_start;
            }
            freq_curr_start = freq_cnt;
            for (i = 0; i < nploty; ++i) {
              for (j = 0; j < nplotx; ++j) {
                if (freq_cnt < numfreq) {
                  polar(k*90, clocks_array[k][freq_cnt],
                        2*j*r + r, 2*i*r + r, r,
                        pdmax[k][freq_cnt]);
                }
                freq_cnt = freq_cnt + 1;
                // cout << "increased f: " << freq_cnt << endl;
              }
            }
            newdraw = false;

          } else if(e.type==ButtonPress) {
            // cout << "Buttonpress found" << endl;
            XClearWindow(dpy, w);
            newdraw = true;
            break;
          }

          else if (e.type==KeyPress) {
            // cout << "keyPress found" << endl;
            // cout << e.xkey.state << endl;
            // cout << XKeycodeToKeysym(dpy,e.xkey.keycode,0) << endl;
            if (XLookupKeysym(&e.xkey, 0) != XK_Alt_L) {
              XClearWindow(dpy, w);
              newdraw = true;
              break;
            }
          }
        } // redrawing loop
      }
    }
  }
} /*manyplots*/
#endif

/*****************************************************************************/
/*
procedure dependence;
var gmin, gmax, fbmin, fbmax, pbdx, x, y: single;
    i, xmax, ymax, xob, yob: word;
{ 
  if (initgr
    ) {
       
         gmin = pdmax[1][1];
         gmax = pdmax[1][1];
         for i = 2 to numfreq do
             {
               if (pdmax[1][i]<gmin
                 ) {
                    gmin = pdmax[1][i];
               if (pdmax[1][i]>gmax
                 ) {
                    gmax = pdmax[1][i]
             };
         xmax = getmaxx;
         ymax = getmaxy;
         x = xmax/8;
         pbdx = deltafreq*3/4*xmax/(freqh-freqd);
         for i = 1 to numfreq do
             {
               y = 5/6*ymax-(pdmax[1][i]-gmin)/(gmax-gmin)
                                                  *2/3*ymax;
               circle(round(x),round(y),3);
               x = x+pbdx
             };
         readkey;
         closegraph
       }
} *//* dependence */

/*****************************************************************************/
int main (int argc, char* argv[]) {

  if (argc == 1) {
    cout << "Error: Missing command line argument (input file name) !" << endl;
    exit(1);
  }

  string arg = argv[1];
  cout << "agv1" << arg << endl;
  string look_for = "nec";
  if (arg.find(look_for) != string::npos) {
    inputNEC(argc, argv);
  } else {
    input(argc, argv);
  }
  if (!scaling) {
    // new(pfc);
    plasyd();
    cgapex();
    for (unsigned int i = 0; i < ndns1; ++i) {
      rightz[i] = right_first[i];
    }
    cgauss();
    // dispose(pfc);
    camplit();
    inimp();
    pabs();
    // if (gd) {graphdistrib()};
    // new(h1);
    get_clock(argv, 0);
    // dispose(h1);
  }
  else {

    deltafreq = (freqh-freqd)/(numfreq-1);
    // new(ph);
    for (unsigned int pordip = 0; pordip < ndip; ++pordip) {
      dhl[pordip] = hl[pordip]*deltafreq;
      dposy[pordip] = posy[pordip]*deltafreq;
      dposz[pordip] = posz[pordip]*deltafreq;
      dar[pordip] = ar[pordip]*deltafreq;
      dhgap[pordip] = hgap[pordip]*deltafreq;
      hl[pordip] = hl[pordip]*freqd;
      posy[pordip] = posy[pordip]*freqd;
      posz[pordip] = posz[pordip]*freqd;
      ar[pordip] = ar[pordip]*freqd;
      hgap[pordip] = hgap[pordip]*freqd;
    };

    for (freq_cnt = 0; freq_cnt < numfreq; ++freq_cnt) {
      cout << fixed << endl;
      cout << "Scaling  " << freqd + freq_cnt*deltafreq << endl;
      if (lambda != 1) {
        cout << " f=" << (freqd + freq_cnt*deltafreq)*299800/lambda
             << " MHz" << endl;
      }
      cout << "ndip: " << ndip << "  nsub: " << nsub;
      cout << "  ngaux: " << ngaux << "  ngaufi: " << ngaufi << endl;
      cout << " num  hl     posy    posz      ar        hgap      v.re   v.im";
      cout << endl;
      for (unsigned int i = 0; i < ndip; i++) {
        cout << setw(3) << i + 1;
        cout << fixed;
        cout.precision(4);
        cout << " " << setw(3) << hl[i];
        cout.precision(4);
        cout << " " << setw(7) << posy[i] << " " <<  setw(7) << posz[i];
        cout << scientific;
        cout << " " << setw(7) << ar[i] << " " << hgap[i];
        cout << fixed;
        cout.precision(2);
        cout << " " << setw(6) << v[i].real();
        cout << " " << setw(6) << v[i].imag() << endl;
      }
      cout << scientific;

      // new(pfc);
      plasyd();
      cgapex();
      for (unsigned int i = 0; i < ndns1; ++i) {
        rightz[i] = right_first[i];
      }
      cgauss();
      // dispose(pfc);
      camplit();
      inimp();
      pabs();
      // if (gd) {graphdistrib();}
      get_clock(argv, freq_cnt);
      for (unsigned int pordip = 0; pordip < ndip; ++pordip) {
        hl[pordip] = hl[pordip]+dhl[pordip];
        posy[pordip] = posy[pordip]+dposy[pordip];
        posz[pordip] = posz[pordip]+dposz[pordip];
        ar[pordip] = ar[pordip]+dar[pordip];
        hgap[pordip] = hgap[pordip]+dhgap[pordip];
      }
    }

    cout << endl;
    cout << "               E-plane            H-plane        Input impedance"
         << endl;
    cout << "Scaling   Thetamax  Gainmax  Thetamax Gainmax      Re   Im   "
         << endl;
    for (unsigned int freq_cnt = 0; freq_cnt < numfreq; ++freq_cnt) {
      cout << fixed << setprecision(2) << setw(6);
      cout << freqd+freq_cnt*deltafreq;
      cout << " " << setw(9) << pthetamax[0][freq_cnt];
      cout << " " << setw(9) << setprecision(3) << 10*log(pdmax[0][freq_cnt])/log(10);
      cout << " " << setw(8) << setprecision(2) << pthetamax[1][freq_cnt];
      cout << " " << setw(9) << setprecision(3) << 10*log(pdmax[1][freq_cnt])/log(10);
      cout << " " << setw(8) << pinimped[freq_cnt].real();
      cout << " " << setw(8) << pinimped[freq_cnt].imag() << endl;
    };
#ifdef _WIN32
#else
    manyplots();
#endif
    // dependence();
    // dispose(ph);
  }
}
