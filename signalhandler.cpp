#include "signalhandler.h"
#include "exceptions.h"
#include <sstream>

namespace ProcessControllerSignalHandler
{

void pushSignalHandler(bool save_handlers)
{
  SignalHandlersStorage dummy_storage;
  SignalHandlersStorage& storage_ref = save_handlers ? handlers_storage : dummy_storage;

  storage_ref.handler_SIGINT = signal(SIGINT, processControllerSignalHandler);
  storage_ref.handler_SIGILL = signal(SIGILL, processControllerSignalHandler);
  storage_ref.handler_SIGFPE = signal(SIGFPE, processControllerSignalHandler);
  storage_ref.handler_SIGSEGV = signal(SIGSEGV, processControllerSignalHandler);
  storage_ref.handler_SIGTERM = signal(SIGTERM, processControllerSignalHandler);
  storage_ref.handler_SIGABRT = signal(SIGABRT, processControllerSignalHandler);
}

void popSignalHandler()
{
  signal(SIGINT, handlers_storage.handler_SIGINT);
  signal(SIGILL, handlers_storage.handler_SIGILL);
  signal(SIGFPE, handlers_storage.handler_SIGFPE);
  signal(SIGSEGV, handlers_storage.handler_SIGSEGV);
  signal(SIGTERM, handlers_storage.handler_SIGTERM);
  signal(SIGABRT, handlers_storage.handler_SIGABRT);

}

void processControllerSignalHandler(int sig)
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

  std::stringstream sstm;
  sstm << "ERROR: System failure. Rank: " << error_handler_MPI_rank 
       << ", Signal: " << sig << " (" << msgsig << ")";
  throw ProcessTaskException(sstm.str());
}

} // namespace ProcessControllerSignalHandler
