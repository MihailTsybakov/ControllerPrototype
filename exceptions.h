#ifndef EXCEPTIONS
#define EXCEPTIONS

#include <exception>
#include <string>

#define DECLARE_EXCEPTION(ClassName, BaseClassName) \
class ClassName: public BaseClassName { \
private: \
  int err_code; \
public: \
  ClassName(const char* msg): BaseClassName(msg) {} \
  ClassName(const char* msg, int err_code): BaseClassName(msg), err_code(err_code) {} \
  ClassName(const std::string& msg): BaseClassName(msg.c_str()) {} \
  ClassName(const std::string& msg, int err_code): BaseClassName(msg.c_str()), err_code(err_code) {} \
   \
  int errCode() const noexcept { return err_code; } \
};

DECLARE_EXCEPTION(ProcessControllerException, std::runtime_error)
DECLARE_EXCEPTION(ProcessTaskException, std::runtime_error)

#endif//EXCEPTIONS