#ifndef PROCESS_CONTROLLER_SIGNAL_HANDLER
#define PROCESS_CONTROLLER_SIGNAL_HANDLER

#include <csignal>
#include <iostream>
#include <stack>

static int error_handler_MPI_rank = 0;

void ProcessControllerSignalHandler(int sig);
void PushSignalHandler(bool save_handlers = true);
void PopSignalHandler();

typedef void(*SignalHandlerPointer)(int);

typedef struct SignalsHandlerPointer {
	SignalHandlerPointer Handler_SIGINT;
	SignalHandlerPointer Handler_SIGILL;
	SignalHandlerPointer Handler_SIGFPE;
	SignalHandlerPointer Handler_SIGSEGV;
	SignalHandlerPointer Handler_SIGTERM;
	SignalHandlerPointer Handler_SIGABRT;
	SignalHandlerPointer Handler_SIGABRT_COMPAT;
	SignalsHandlerPointer() {
		Handler_SIGINT = nullptr;
		Handler_SIGILL = nullptr;
		Handler_SIGFPE = nullptr;
		Handler_SIGSEGV = nullptr;
		Handler_SIGTERM = nullptr;
		Handler_SIGABRT = nullptr;
		Handler_SIGABRT_COMPAT = nullptr;
	}
}SignalsHandlerPointer;

static SignalsHandlerPointer handlers;

#endif//PROCESS_CONTROLLER_SIGNAL_HANDLER
