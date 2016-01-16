#include	<signal.h>
#include	<sys/types.h>
#include	<unistd.h>
#include	<sys/wait.h>
#include        <errno.h>

#define USRSHELL	"/usr/bin/ksh"
#define SINGLESHELL	"/sbin/sh"

int 	
system(const char *string)
{
	sigset_t old,new;
	int     status, w;
	pid_t   pid;
	void	(*oldint)();
	void	(*oldquit)();
	char   *shell, *argstr;
	extern  char _DYNAMIC[];

	if(string == NULL) {
	  if (access(USRSHELL, X_OK) && access(SINGLESHELL, X_OK))
	      return(0);
	  else
	      return(1);
	}

	sigemptyset( &new );
	sigemptyset( &old );
	sigaddset( &new, SIGCHLD);
	sigprocmask( SIG_BLOCK, &new, &old);
	oldint = signal(SIGINT, SIG_IGN);
	oldquit = signal(SIGQUIT, SIG_IGN);

	if ((pid = vfork()) == -1) {
		/* vfork() will have set errno */
		sigprocmask( SIG_SETMASK, &old, NULL);
		signal(SIGINT, oldint);
		signal(SIGQUIT, oldquit);
		return((string) ? -1 : 0);
	}
	if (!pid) { /* child */
	        sigprocmask( SIG_SETMASK, &old, NULL);
		signal(SIGINT, oldint);
		signal(SIGQUIT, oldquit);
		shell = (access(USRSHELL, X_OK)==0 && _DYNAMIC)
				? USRSHELL : SINGLESHELL;

		argstr = (string) ? (char *) string : "";

		(void) execl(shell, "sh", "-c", argstr, 0);

		_exit(127);
	}

	do {
	    w = waitpid( pid, &status, 0);
	} while (w == -1 && errno == EINTR);

	sigprocmask( SIG_SETMASK, &old, NULL);
	signal(SIGINT, oldint);
	signal(SIGQUIT, oldquit);

	if (w == -1)
		return((string) ? -1 : 0);

	return((string) ? status : (status == 0));
}

