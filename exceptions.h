#ifndef EXCEPTIONS
#define EXCEPTIONS

#include <exception>
#include <string>

class ProcessControllerException : public std::exception
{
private:
	std::string error_msg;
	int error_code;
public:
	ProcessControllerException(std::string error_message, int err_code) : error_msg(error_message), error_code(err_code) {}
	int errCode() const noexcept { return error_code; }
	const char* what() const noexcept override { return error_msg.c_str(); }
};

class ProcessTaskException : public std::exception
{
private:
	std::string error_msg;
	int error_code;
public:
	ProcessTaskException(std::string error_message, int err_code) : error_msg(error_message), error_code(err_code){}
	int errCode() const noexcept { return error_code; }
	const char* what() const noexcept override { return error_msg.c_str(); }
};

#endif//EXCEPTIONS