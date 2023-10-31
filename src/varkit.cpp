#include <zlib.h>
#include <string>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

#define LENGTH 51200

std::string next_chunk(gzFile file) {
  int bytes_read;
  char buffer[LENGTH];

  bytes_read = gzread(file, buffer, LENGTH - 1);
  buffer[bytes_read] = '\0';

  return reinterpret_cast<char*>(buffer);
}

void str_split(std::string& str, char& split, std::vector<std::string>& result){
  unsigned int start = 0;
  unsigned int i = 0;

  for(i = 1; i < str.size(); i++){
    if(str[i] == split){
      std::string temp = str.substr(start, i - start);
      result.push_back(temp);
      start = i + 1;
    }
  }

  // Handle last
  std::string temp = str.substr(start, i - start);
  result.push_back(temp);
}

// [[Rcpp::export(name=".str_splitC")]]
CharacterVector str_split(std::string str, char sep) {
  std::vector<std::string> lines;
  str_split(str, sep, lines);
  return wrap(lines);
}

std::string substr_till(std::string& str, char& till, unsigned int& n) {
  unsigned int count = 0;
  unsigned int i;
  unsigned int end = str.length();
  for (i = 0; i < end; i++) {
    if (str[i] != till) {
      continue;
    }
    count++;
    if (count >= n) {
      break;
    }
  }
  return str.substr(0, i);
}

std::string get_chrom(std::string& str) {
  if (str[1] == '#' && str[0] == '#') {
    return "##";
  }
  char till = '\t';
  unsigned int count = 1;
  return substr_till(str, till, count);
}

// [[Rcpp::export(name=".index_vcfC")]]
DataFrame index_vcf(std::string x) {
  gzFile file = gzopen(x.c_str(), "r");

  if (! file) {
    stop("Could not open '" + x + "'! Error: " + strerror(errno));
  }

  // Vectors to store results
  std::vector<std::string> id;
  std::vector<int> n_lines;
  std::string last_line = "";
  int err;

  // Looping variables
  std::string current = "";
  unsigned int n_current = 0;

  while (1) {
    checkUserInterrupt();
    std::string str = next_chunk(file);
    str = last_line + str;

    char split = '\n';
    std::vector<std::string> lines;
    str_split(str, split, lines);

    // Skip the last line as it might be incomplete
    int i_stop = lines.size() - 1;
    for (int i = 0; i < i_stop; i++) {
      std::string chrom = get_chrom(lines[i]);
      // Rcout << chrom << std::endl;

      if (chrom == current) {
        n_current++;
        continue;
      }

      if (n_current > 0) {
        id.push_back(current);
        n_lines.push_back(n_current);
      }

      current = chrom;
      n_current = 1;
    }

    last_line = lines[lines.size() - 1];

    if (gzeof(file)) {
      break;
    } else {
      const char * error_string = gzerror(file, & err);
      if (err) {
        stop(std::string("Error: ") + error_string);
      }
    }
  }

  // Handle last
  id.push_back(current);
  n_lines.push_back(n_current);

  if (last_line.length() != 0) {
    warning("The file is not terminated by a new line, the index might be off.");
  }

  gzclose(file);

  return DataFrame::create(Named("id") = id, Named("n_lines") = n_lines);
}
