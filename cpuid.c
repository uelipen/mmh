#include <cpuset.h>
int cpuid_()
{
	return cpu_get_current();
}

#include <numa.h>

int cpu_get_rad_()
{
	return cpu_get_rad(cpu_get_current());
}
