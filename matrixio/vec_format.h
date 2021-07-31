#ifndef VEC_FORMAT_H
#define VEC_FORMAT_H

namespace MATRIXIO
{
template<class T>
void read_vector(std::ifstream& file, bool binary, MKL_INT full_size, std::vector<T>& result, bool vecInt32 = false, int first = 0, int last = -1)
{
  if (last == -1)
    last = static_cast<int>(full_size - 1);
  int size = last - first + 1;

  result.resize(size);
  if (binary) {
    if (vecInt32) {
      int temp = static_cast<int>(full_size);
      boost::scoped_array<T> array(new T[full_size]);
      file.read(reinterpret_cast<char*>(array.get()), sizeof(T) * full_size);
      for (MKL_INT i = 0; i < size; ++i)
        result[i] = array[first + i];
    }
    else {
      boost::scoped_array<T> array(new T[full_size]);
      file.read(reinterpret_cast<char*>(array.get()), sizeof(T) * full_size);
      for (MKL_INT i = 0; i < size; ++i)
        result[i] = array[first + i];
    }

  }
  else {
    T tmp;
    for (MKL_INT i = 0; i < full_size; ++i) {
      if (i >= first && i <= last)
        file >> result[i - first];
      else
        file >> tmp;
    }
  }
}

template<class T>
void read_vector_with_shift(std::ifstream& file, bool binary, MKL_INT full_size, std::vector<T>& result, int& shift, bool vecInt32 = false, int first = 0, int last = -1)
{
  if (last == -1)
    last = full_size - 1;
  int size = last - first + 1;

  result.resize(size);
  if (binary) {
    if (vecInt32) {
      int temp = full_size;
      boost::scoped_array<T> array(new T[full_size]);
      file.read(reinterpret_cast<char*>(array.get()), sizeof(T) * full_size);
      shift = array[first] - array[0];
      for (MKL_INT i = 0; i < size; ++i)
        result[i] = array[first + i] - shift;
    }
    else {
      boost::scoped_array<T> array(new T[full_size]);
      file.read(reinterpret_cast<char*>(array.get()), sizeof(T) * full_size);
      shift = array[first] - array[0];
      for (MKL_INT i = 0; i < size; ++i)
        result[i] = array[first + i] - shift;
    }

  }
  else {
    T tmp;
    shift = 0;
    for (MKL_INT i = 0; i < full_size; ++i) {
      file >> tmp;
      if (i == 0)
        shift = tmp;
      if (i == first)
        shift = tmp - shift;
      if (i >= first && i <= last)
        result[i - first] = tmp - shift;
    }
  }
}

template<class T>
void read_vector_file(std::string& file_name,
  std::vector<T>& result,
  MKL_INT& vector_size,
  bool binary, bool vecInt32 = false)
{
  std::ifstream file(file_name.c_str(), std::ifstream::in | std::ifstream::binary);
  if (!file.is_open())
  {
    std::cerr << "   ERROR: Cannot open file " << file_name << ".\n";
    return;
  }

  if (vecInt32) {
    int size;
    if (binary)
      file.read(reinterpret_cast<char*>(&size), sizeof(size));
    else
      file >> size;
    vector_size = size;
    read_vector(file, binary, size, result, vecInt32);
  }
  else {
    MKL_INT size;
    if (binary)
      file.read(reinterpret_cast<char*>(&size), sizeof(size));
    else
      file >> size;
    vector_size = size;
    read_vector(file, binary, size, result, vecInt32);
  }

}

template<class T>
void save_vector(std::ofstream& file, bool binary, const std::vector<T>& array)
{
  if (binary)
  {
    file.write(reinterpret_cast<const char*>(&array[0]), sizeof(T) * array.size());
  }
  else
  {
    file << std::setw(20) << std::setprecision(20);
    for (unsigned i = 0; i < array.size(); ++i)
    {
      file << std::scientific << array[i] << "\n";
    }
  }
}

template<class T>
void save_vector_file(std::ofstream& file, bool binary, const std::vector<T>& array)
{
  if (binary)
  {
    const MKL_INT size = (MKL_INT)array.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
  }
  else
  {
    file << array.size() << "\n";
  }
  save_vector(file, binary, array);
  file.flush();
}

}

#endif 
