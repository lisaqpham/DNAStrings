//DNA Strings - Data Structures Assignment 1
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

// declare variables
int lines;
int sum;
double mean;
double diffSum;
double sqDiff;
double variance;
double stdev;
double a=0, c=0, t=0, g=0;
double aa=0, ac=0, at=0, ag=0;
double ca=0, cc=0, ct=0, cg=0;
double ta=0, tc=0, tt=0, tg=0;
double ga=0, gc=0, gt=0, gg=0;
char nucleotide;
char bigram;
string dnaString;
string newDNA;
string answer;
string newFile;
bool repeat = true;

// resets file "cursor"
void resetCursor(istream& file){
  file.clear();
  file.seekg(0);
}

// generates lengths for new DNA strings based on Gaussian Distribution
int generateLength() {
  double x, y, C, D;
  x = (double)rand() / RAND_MAX;
  y = (double)rand() / RAND_MAX;

  C = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
  D = sqrt(variance) * C + mean;

  return (int)round(D);
}

// generates nucleotides for the new DNA strings based on prior probability
char generateNucleo() {
  double r = (double)rand() / RAND_MAX;
  if (r <= a) {
    return 'a';
  } else if (r <= (a+c)) {
    return 'c';
  } else if (r <= (a+c+t)) {
    return 't';
  } else {
    return 'g';
  }
}

// computes sum
void findSum(istream& file) {
  lines = 0;
  sum = 0;
  while(file >> dnaString) {
    lines++;
    sum += dnaString.length();
  }
}

// computes mean
void findMean(istream& file) {
  mean = 0;
  mean = (double)sum / lines;
}

// computes variance
void findVariance(istream& file) {
  diffSum = 0;
  resetCursor(file);
  while(file >> dnaString) {
    sqDiff = (dnaString.length() - mean);
    sqDiff = sqDiff * sqDiff;
    diffSum += sqDiff;
  }
  variance = diffSum / lines;
}

// computes standard deviation
void findStDev(istream& file){
  stdev = sqrt(variance);
}

// computes relative probability
void relativeProbability(istream& file){
  resetCursor(file);
  while (file >> nucleotide) {
    if (tolower(nucleotide) == 'a') {
      a++;
    } else if (tolower(nucleotide) == 'c') {
      c++;
    } else if (tolower(nucleotide) == 't') {
      t++;
    } else if (tolower(nucleotide) == 'g') {
      g++;
    }
  }
  a /= sum;
  c /= sum;
  t /= sum;
  g /= sum;
}

// computes bigram probability
void bigramProbability(istream& file){
  resetCursor(file);
  while (file >> bigram) {
    if (tolower(bigram) == 'a') {
      if (tolower(file.peek()) == 'a'){
        aa++;
      } else if (tolower(file.peek()) == 'c'){
        ac++;
      } else if (tolower(file.peek()) == 't'){
        at++;
      } else if (tolower(file.peek()) == 'g'){
        ag++;
      }
    } else if (tolower(bigram) == 'c') {
      if (tolower(file.peek()) == 'a'){
        ca++;
      } else if (tolower(file.peek()) == 'c'){
        cc++;
      } else if (tolower(file.peek()) == 't'){
        ct++;
      } else if (tolower(file.peek()) == 'g'){
        cg++;
      }
    } else if (tolower(bigram) == 't') {
      if (tolower(file.peek()) == 'a'){
        ta++;
      } else if (tolower(file.peek()) == 'c'){
        tc++;
      } else if (tolower(file.peek()) == 't'){
        tt++;
      } else if (tolower(file.peek()) == 'g'){
        tg++;
      }
    } else if (tolower(bigram) == 'g') {
      if (tolower(file.peek()) == 'a'){
        ga++;
      } else if (tolower(file.peek()) == 'c'){
        gc++;
      } else if (tolower(file.peek()) == 't'){
        gt++;
      } else if (tolower(file.peek()) == 'g'){
        gg++;
      }
    }
  }
  aa /= (sum - lines);
  ac /= (sum - lines);
  at /= (sum - lines);
  ag /= (sum - lines);

  ca /= (sum - lines);
  cc /= (sum - lines);
  ct /= (sum - lines);
  cg /= (sum - lines);

  ta /= (sum - lines);
  tc /= (sum - lines);
  tt /= (sum - lines);
  tg /= (sum - lines);

  ga /= (sum - lines);
  gc /= (sum - lines);
  gt /= (sum - lines);
  gg /= (sum - lines);
}

// outputs statistics for DNA strings file
void outputStats(ofstream& file) {
  if (file.is_open()) {
    file << "Name: Lisa Pham" << endl;
    file << "ID: 2338933" << endl;
    file << "Assignment: DNA Strings C++ Review" << endl;
    file << "Class: CPSC 350-02" << endl;
    file << "File Name: " << newFile << endl << endl;

    file << "Lines: " << lines << endl;
    file << "Sum: " << sum << endl;
    file << "Mean: " << mean << endl;
    file << "Variance: " << variance << endl;
    file << "StDev: " << stdev << endl << endl;

    file << "Relative Probability: " << endl;
    file << "A: " << a << endl;
    file << "C: " << c << endl;
    file << "T: " << t << endl;
    file << "G: " << g << endl << endl;

    file << "Bigram Probability: " << endl;
    file << "AA: " << aa << endl;
    file << "AC: " << ac << endl;
    file << "AT: " << at << endl;
    file << "AG: " << ag << endl << endl;

    file << "CA: " << ca << endl;
    file << "CC: " << cc << endl;
    file << "CT: " << ct << endl;
    file << "CG: " << cg << endl << endl;

    file << "TA: " << ta << endl;
    file << "TC: " << tc << endl;
    file << "TT: " << tt << endl;
    file << "TG: " << tg << endl << endl;

    file << "GA: " << ga << endl;
    file << "GC: " << gc << endl;
    file << "GT: " << gt << endl;
    file << "GG: " << gg << endl << endl;
  } else {
    cout << "Error: Cannot write to file." << endl;
  }
}

void generate1000(ofstream& file){
  file << "Generated Gaussian DNA (1000 strings)" << endl;
  for(int i = 0; i < 1000; i++) {
    newDNA = "";
    for(int j = 0; j < generateLength(); j++) {
      newDNA += generateNucleo();
    }
    file << newDNA << endl;
  }
  file << endl;
}


// Main Method
int main(int argc, char** argv) {

    ofstream outfile("lisapham.out");
    cout << "Please type your file name: " << endl;
    cin >> newFile;

    // loop for multiple file inputs
    while(repeat) {
      // file open for reading / writing
      ifstream infile(newFile);

      // check if file exists
      if (!infile) {
        cout << "Error: Cannot read file." << endl;
        cout << "Please re-type the name of your file." << endl;
        cin >> newFile;
        ifstream infile(newFile);
        resetCursor(infile);
      }

      // basic stats
      findSum(infile);
      findMean(infile);
      findVariance(infile);
      findStDev(infile);

      // probability stats
      relativeProbability(infile);
      bigramProbability(infile);

      // ouput statistics to file
      outputStats(outfile);

      // generating DNA Strings
      generate1000(outfile);

      // resets infile
      infile.close();
      infile.clear();


      // prompt for new file
      cout << "Would you like to process another list? (y/n)" << endl;
      cin >> answer;

      if (answer == "y" || answer == "yes") {
        cout << "Please type the name of your file: " << endl;
        cin >> newFile;
      } else {
        cout << "Thank you for using this program." << endl;
        // close and exit
        outfile.close();
        repeat = false;
        return 0;
      }
    }
}
