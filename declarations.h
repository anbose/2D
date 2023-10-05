void send_boundary(const int N, const int count, const int bnum, const double *localrho, double *leftboundary, double *rightboundary);
void send_boundary_from_global(const int N, const int bnum, const int *firstindex, const double *globalrho, double *leftboundary, double *rightboundary);
