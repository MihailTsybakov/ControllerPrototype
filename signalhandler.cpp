#include "signalhandler.h"

void PushSignalHandler(bool save_handlers)
{
	if (save_handlers == true)
	{
		handlers.Handler_SIGINT = signal(SIGINT, ProcessControllerSignalHandler);
		handlers.Handler_SIGILL = signal(SIGILL, ProcessControllerSignalHandler);
		handlers.Handler_SIGFPE = signal(SIGFPE, ProcessControllerSignalHandler);
		handlers.Handler_SIGSEGV = signal(SIGSEGV, ProcessControllerSignalHandler);
		handlers.Handler_SIGTERM = signal(SIGTERM, ProcessControllerSignalHandler);
		handlers.Handler_SIGABRT = signal(SIGABRT, ProcessControllerSignalHandler);
	}
	else
	{
		signal(SIGINT, ProcessControllerSignalHandler);
		signal(SIGILL, ProcessControllerSignalHandler);
		signal(SIGFPE, ProcessControllerSignalHandler);
		signal(SIGSEGV, ProcessControllerSignalHandler);
		signal(SIGTERM, ProcessControllerSignalHandler);
		signal(SIGABRT, ProcessControllerSignalHandler);
	}
}

void PopSignalHandler()
{
	handlers.Handler_SIGINT = signal(SIGINT, handlers.Handler_SIGINT);
	handlers.Handler_SIGILL = signal(SIGILL, handlers.Handler_SIGILL);
	handlers.Handler_SIGFPE = signal(SIGFPE, handlers.Handler_SIGFPE);
	handlers.Handler_SIGSEGV = signal(SIGSEGV, handlers.Handler_SIGSEGV);
	handlers.Handler_SIGTERM = signal(SIGTERM, handlers.Handler_SIGTERM);
	handlers.Handler_SIGABRT = signal(SIGABRT, handlers.Handler_SIGABRT);

}

void ProcessControllerSignalHandler(int sig)
{
	const char* msgsig = "Signal unknown";
	switch (sig)
	{
	case 1:
		msgsig = "Hangup";
		break;
	case 2:
		msgsig = "Terminal interrupt signal (CTRL-C)";
		break;
	case 3:
		msgsig = "Terminal quit signal";
		break;
	case 4:
		msgsig = "Illegal Instruction Processor";
		break;
	case 5:
		msgsig = "Trace/breakpoint trap";
		break;
	case 6:
		msgsig = "The process is complete function is called Abort() ";
		break;
	case 14:
		msgsig = "Alarm clock";
		break;
	case 10:
		msgsig = "Access to an undefined portion of a memory object";
		break;
	case 18:
		msgsig = "Child process terminated, stopped";
		break;
	case 25:
		msgsig = "Continue executing, if stopped";
		break;
	case 8:
		msgsig = "Erroneous arithmetic operation";
		break;
	case 9:
		msgsig = "Kill (cannot be caught or ignored)";
		break;
	case 13:
		msgsig = "Write on a pipe with no one to read it";
		break;
	case 11:
		msgsig = "Invalid memory reference";
		break;
	case 15:
		msgsig = "Termination signal";
		break;
	case 20:
		msgsig = "Terminal stop signal";
		break;
	case 26:
		msgsig = "Background process attempting read";
		break;
	case 27:
		msgsig = "Background process attempting write";
		break;
	case 16:
		msgsig = "User-defined signal 1";
		break;
	case 17:
		msgsig = "User-defined signal 2";
		break;
	case 22:
		msgsig = "An unhandled exception was not in the main thread ";
		break;
	case 29:
		msgsig = "Profiling timer expired";
		break;
	case 12:
		msgsig = "Bad system call";
		break;
	case 21:
		msgsig = "High bandwidth data is available at a socket";
		break;
	case 28:
		msgsig = "Virtual timer expired";
		break;
	case 30:
		msgsig = "CPU time limit exceeded";
		break;
	case 31:
		msgsig = "File size limit exceeded";
		break;
	default:
		break;
	}

	std::cout << "ERROR: System failure. Rank: " << error_handler_MPI_rank << ", Signal: " << sig << " (" << msgsig << ")" << std::endl;
}