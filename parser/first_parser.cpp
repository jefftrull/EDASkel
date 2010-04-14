#include <iostream>
#include <string>
#include <vector>

#include "defparser.h"

using namespace std;

int main(int argc, char **argv) {

  std::string filecontents("DESIGN test ;\nVERSION 1.211 ;\nEND DESIGN\n");
  std::string::const_iterator beg = filecontents.begin();
  std::string::const_iterator end = filecontents.end();

  defparser<string::const_iterator> defp;
  bool pass = phrase_parse(beg, end, defp, boost::spirit::qi::space);

  if (!pass)
    cerr << "failed!\n";
  else
    cerr << "passed!\n";

  return 0;
}
