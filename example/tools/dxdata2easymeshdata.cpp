/**
 * 将 OpenDX 格式的数据文件转换为 EasyMesh 格式的数据文件
 * 
 */

#include "Geometry.h"
#include "EasyMesh.h"
#include "TemplateElement.h"

#define DIM 2
#define DOW 2

using namespace AFEPack;

SimplestMesh<DIM,DOW> simp_mesh;
EasyMesh easy_mesh;
std::vector<double> data;

void readPoint(std::istream&);
void readConnection(std::istream&);
void readData(std::istream&);

int main(int argc, char * argv[])
{
  if (argc != 3 && argc != 4) {
    std::cout << "Usage: " << argv[0]
	      << " input_file output_file have_data\n"
	      << "\tinput_file: Input file with Open DX data file format\n"
	      << "\toutput_file: Output file with EasyMesh data file format\n"
	      << "\thave_data: if the file have data (!=0) or have only mesh info (==0)\n"
	      << std::endl;
    return 1;
  }
  std::string infile(argv[1]);
  std::string outfile(argv[2]);

  std::vector<TemplateGeometry<DIM> > template_geometry(1);
  template_geometry[0].readData("triangle.tmp_geo");
  simp_mesh.setTemplateGeometry(template_geometry);

  std::ifstream is(infile.c_str());
  readPoint(is);
  readConnection(is);

  simp_mesh.generateMesh(easy_mesh);
  easy_mesh.writeData(outfile.c_str());

  if (argc == 3 || atoi(argv[3])) {
    readData(is);
    std::ofstream os((outfile + ".dat").c_str());
    for (int i = 0;i < data.size();i ++) {
      os << data[i] << "\n";
    }
    os.close();
  }

  return 0;
}

void readPoint(std::istream& is)
{
  std::cout << "\tReading in points data ... " << std::flush;
  int n_point;
  char buffer[256];
  std::string word;
  do {
    is.get(buffer, 256, ' ');
    word = buffer;
    if (word == "item") {
      is >> n_point;
      break;
    }
    is.get(buffer[0]);
  } while (1);
  is.get(buffer, 256);
  std::cout << n_point << " total, ... " << std::flush;
  simp_mesh.point().resize(n_point);
  for (int i = 0;i < n_point;i ++) {
    is >> simp_mesh.point(i);
  }
  std::cout << "OK!" << std::endl;
};

void readConnection(std::istream& is)
{
  std::cout << "\tReading in connections data ... " << std::flush;
  int n_element;
  char buffer[256];
  std::string word;
  do {
    is.get(buffer, 256, ' ');
    word = buffer;
    if (word == "item") {
      is >> n_element;
      break;
    }
    is.get(buffer[0]);
  } while (1);
  is.get(buffer, 256);
  simp_mesh.element().resize(n_element);
  std::cout << n_element << " total, ... " << std::flush;
  for (int i = 0;i < n_element;i ++) {
    simp_mesh.element(i).template_element = 0;
    simp_mesh.elementVertex(i).resize(DIM + 1);
    for (int j = 0;j < DIM + 1;j ++) {
      is >> simp_mesh.elementVertex(i, j);
    }
  }
  std::cout << "OK!" << std::endl;
}

void readData(std::istream& is)
{
  std::cout << "\tReading in funciton data ... " << std::flush;
  int n_data;
  char buffer[256];
  std::string word;
  do {
    is.get(buffer, 256, ' ');
    word = buffer;
    if (word == "item") {
      is >> n_data;
      break;
    }
    is.get(buffer[0]);
  } while (1);
  is.get(buffer, 256);
  std::cout << n_data << " total, ... " << std::flush;
  Assert(n_data == simp_mesh.n_point(), ExcInternalError());
  data.resize(n_data);
  for (int i = 0;i < n_data;i ++) {
    is >> data[i];
  }
  std::cout << "OK!" << std::endl;
}

/**
 * end of file
 * 
 */
