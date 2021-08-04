#ifndef PROCESS_CONTROLLER_SIGNAL_HANDLER
#define PROCESS_CONTROLLER_SIGNAL_HANDLER

#include <csignal>
#include <iostream>
#include <stack>

namespace ProcessControllerSignalHandler
{

static int error_handler_MPI_rank = 0;

void processControllerSignalHandler(int sig);
void pushSignalHandler(bool save_handlers = true);
void popSignalHandler();

typedef void(*SignalHandlerPtr)(int);

typedef struct SignalHandlersStorage {
  SignalHandlerPtr handler_SIGINT;
  SignalHandlerPtr handler_SIGILL;
  SignalHandlerPtr handler_SIGFPE;
  SignalHandlerPtr handler_SIGSEGV;
  SignalHandlerPtr handler_SIGTERM;
  SignalHandlerPtr handler_SIGABRT;
  SignalHandlerPtr handler_SIGABRT_COMPAT;
  SignalHandlersStorage() {
    handler_SIGINT = nullptr;
    handler_SIGILL = nullptr;
    handler_SIGFPE = nullptr;
    handler_SIGSEGV = nullptr;
    handler_SIGTERM = nullptr;
    handler_SIGABRT = nullptr;
    handler_SIGABRT_COMPAT = nullptr;
  }
} SignalHandlersStorage;

static SignalHandlersStorage handlers_storage;

} // namespace ProcessControllerSignalHandler

#endif//PROCESS_CONTROLLER_SIGNAL_HANDLER
