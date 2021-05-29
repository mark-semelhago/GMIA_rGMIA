classdef pardiso < handle
    % PARDISO class implementing the interface ParDiSo-MATLAB.
    %
    %   This class must be intended to be a "static" class, i.e., it is
    %   implemented as a handle. This means that only one single pardiso object
    %   can be created in each MATLAB instance.
    %
    %
    %   The class contains the following properties:
    %   Public:
    %       iparm  - 64-by-1 array of ParDiSo internal integer parameters
    %       dparm  - 64-by-1 array of ParDiSo internal double parameters
    %       solver -  1-by-1 scalar selecting ParDiSo direct/iterative solver
    %
    %
    %   The class contains the following methods:
    %   Public:
    %       (constructor) - constructor of the class. Allocates different
    %                       instances of ParDiSo
    %       factorize     - factorizes a given matrix
    %       solve         - solves a linear system with many right-hand sides
    %       invert        - computes the selected inverse of a given matrix
    %       
    %       getiparm      - returns the current array iparm
    %       getdparm      - returns the current array dparm
    %       printiparm    - prints the current array iparm
    %       printdparm    - prints the current array dparm
    %
    %       getcounter    - returns the number of pardiso instances allocated
    %
    %  
    %   Static:
    %       isstructsymm    - checks if a given matrix is structurally symmetric
    %       makestructsymm  - makes a given matrix structurally symmetric
    %       getpardisomtype - returns the ParDiSo mtype of a given matrix
    %
    %
    %
    %   Fabio VERBOSIO, Drosos KOUROUNIS, Olaf SCHENK
    %   Universita` della Svizzera Italiana, October 5th 2017
    %
    %   See also PARDISO.PARDISO.
    %
    
    properties
        iparm
        msglvl
        dparm
    end
    
    properties (Access = private)
        maxfact = 100;
        prds_count
    end

    methods (Access = public)

        function this = pardiso(prdsid, mtype, varargin)
        % PARDISO   constructor for the class ParDiSo
        % 
        %   PARDISO(PRDSID, MTYPE) creates a ParDiSo object for matrix of type
        %   MTYPE. PRDSID indicates the ParDiSo ID, i.e., the index of the
        %   ParDiSo instance that will be allocated. The user is responsible for
        %   handling such ID.
        %
        %   PARDISO(PRDSID, MTYPE, SOLVER) creates a ParDiSo object and forces
        %   ParDiSo to use the indicated solver. If SOLVER = 0, then a sparse
        %   direct solver will be used. If SOLVER = 1, a multi-recursive 
        %   iterative solver will be used instead. The default is 0.
        %
        %   MTYPE indicates the matrix type and must be one of the following
        %   values:
        %
        %       +========+=========================================+
        %       |  MTYPE |        Matrix type                      |
        %       +=====+============================================+
        %       |    1   | real and structurally symmetric         |
        %       |    2   | real and symmetric positive defite      |
        %       |   -2   | and symmetric indefinite                |
        %       +--------+-----------------------------------------+
        %       |    3   | complex and structurally symmetric      |
        %       |    4   | complex and Hermitian positive definite |
        %       |   -4   | complex and Hermitian indefinite        |
        %       |    6   | complex and symmetric                   |
        %       +--------+-----------------------------------------+
        %       |   11   | real and nonsymmetric                   |
        %       |   13   | complex and nonsymmetric                |
        %       +========+=========================================+
        %
        %   Notice that for unstructured matrices (MTYPE = 11, 13) the matrix
        %   type will be automatically changed to structurally symmetric 
        %   (respectively, MTYPE = 1, 3) and the internal function
        %   MAKESTRUCTSYMM will be called. It is recommended to call such
        %   function in advance and set the matrix type to structurally in order
        %   to improve the performance of the code.
        %
        %   See also PARDISO.MAKESTRUCTSYMM.
        %
     
            p = inputParser;

            if nargin ~= 0 && nargin ~= 2 && nargin ~= 3
                error('PARDISO: incorrect number of input parameters.');
            elseif nargin == 0
                help callpardiso;
                return;
            elseif nargin == 3
                solver = varargin{1};
            else
                solver = 0;
            end            

            defaultSolver = 0;   % sparse direct solver (1 -> multi-recursive iter. solver)

            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'mtype',  @(x) isnumeric(x) && isscalar(x) && (round(x) == x));
            addOptional(p, 'solver', defaultSolver, @(x) x == 0 || x == 1);
            parse(p, prdsid, mtype, solver);

            this.msglvl  = 0;

            ompnt = getenv('OMP_NUM_THREADS');
            if isempty(ompnt)
                error('Set OMP_NUM_THREADS to 1 (or more) thread(s) for parallel execution. Quitting...');
            else
                this.iparm(3) = str2num(ompnt);
            end

            if mtype == 11      % Pardiso cannot treat real general
                mtype = 1;      % => make it real struct-symm.
                fprintf('Warning: ParDiSo will handle the matrix as real, structurally symmetric. Number of nonzeros may increase.\n');
            end
            if mtype == 13      % Pardiso cannot treat complex general
                mtype = 3;      % => make it complex struct-symm.
                fprintf('Warning: ParDiSo will handle the matrix as real, structurally symmetric. Number of nonzeros may increase.\n');
            end

            [this.iparm, this.dparm, err] = callpardiso(prdsid, mtype, solver);
            this.chkerr(err);

            this.iparm(12) = 1;    % Solve systems with the transpose of A;
            this.prds_count = prdsid;
        end

        function this = verbose(this)
            this.msglvl = 1;
        end    

        function this = quiet(this)
            this.msglvl = 0;
        end    

        function factorize(this, prdsid, A, varargin)
        % FACTORIZE factorizes a matrix with PARDISO
        %
        %    FACTORIZE(PRDS, PRDSID, A) computes the factorization of the sparse
        %    matrix A on the ParDiSo instance indicated by PRDSID. PRDS must be 
        %    a pardiso object created calling the constructor of the class.
        %
        %    FACTORIZE(PRDS, PRDSID, A, MNUM) computes the factorization of the
        %    sparse matrix A and matrix number MNUM. It is possible to store up
        %    to MAXFACT matrices sharing the same nonzero pattern on the same
        %    ParDiSo instance.
        %
        %    The factorization must be computed prior to any other operation
        %    such as the solution of a system of the selected inversion. In case
        %    of symmetric matrices (matrix type = 2, -2, 6) only the lower
        %    triangular part of A should be passed as an argument.
        %
        %    See also PARDISO.PARDISO, PARDISO.SOLVE, PARDISO.INVERT.
        %    

            if nargin < 3
                error('Input must be the ParDiSo ID and one sparse matrix.');
            end
            
            mnum = 1;
            if nargin == 4
                mnum = varargin{1};
            end
            
            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y)  && size(A,1) == size(A,2));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, mnum);
            
            err  = 0;
            nrhs = 1;
            b    = [];
            x    = [];
            perm = 0;


            %if issymmetric(A)
            %    A = tril(A);
            %end
            
            phase = 12;

            [~, this.iparm, ~, ~, err, this.dparm] = ...
            callpardiso(prdsid, mnum, phase, ...
                        A, ...
                        perm, nrhs, this.iparm, this.msglvl, ...
                        b, x, err, this.dparm);
            this.chkerr(err);

        end

        function [x] = mldivide(this, prdsid, A, b, varargin)
        % SOLVE solves a linear system with PARDISO
        %
        %    X = SOLVE(PRDS, PRDSID, A, B) solves the linear system A*X = B
        %    using the pardiso instance indicated by PRDSID. PRDS must be a
        %    pardiso object created calling the constructor of the class. The
        %    factorization of the matrix must be computed in advance calling
        %    the function FACTORIZE(PRDS, PRDSID, ...).
        % 
        %    X = SOLVE(PRDS, PRDSID, A, B, MNUM) solves the linear system 
        %    A*X = B for the matrix A relative to matrix number MNUM.
        %
        %    In case of symmetric matrices (matrix type = 2, -2, 6) only the
        %    lower triangular part of A should be passed as an argument.
        %
        %    See also PARDISO.PARDISO, PARDISO.FACTORIZE.
        %

            if nargin < 4
                error('Input must be the ParDiSo ID, one sparse matrix, and the RHS [and the matrix number].');
            end
            
            mnum = 1;
            if nargin == 5
                mnum = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y)  && size(A,1) == size(A,2));
            addRequired(p, 'b',      @(y) isnumeric(y) && size(y,1) == size(A,1));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, b, mnum);

            x = this.solve(prdsid, A, b, mnum);

        end

        function [x] = mrdivide(this, prdsid, A, b, varargin)
        % SOLVE solves a linear system with PARDISO
        %
        %    X = SOLVE(PRDS, PRDSID, A, B) solves the linear system A*X = B
        %    using the pardiso instance indicated by PRDSID. PRDS must be a
        %    pardiso object created calling the constructor of the class. The
        %    factorization of the matrix must be computed in advance calling
        %    the function FACTORIZE(PRDS, PRDSID, ...).
        % 
        %    X = SOLVE(PRDS, PRDSID, A, B, MNUM) solves the linear system 
        %    A*X = B for the matrix A relative to matrix number MNUM.
        %
        %    In case of symmetric matrices (matrix type = 2, -2, 6) only the
        %    lower triangular part of A should be passed as an argument.
        %
        %    See also PARDISO.PARDISO, PARDISO.FACTORIZE.
        %

            if nargin < 4
                error('Input must be the ParDiSo ID, one sparse matrix, and the RHS [and the matrix number].');
            end
            
            mnum = 1;
            if nargin == 5
                mnum = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y)  && size(A,1) == size(A,2));
            addRequired(p, 'b',      @(y) isnumeric(y) && size(y,2) == size(A,1));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, b, mnum);

            this.iparm(12) = 0;
            x = this.solve(prdsid, A, b.', mnum).';
            this.iparm(12) = 1;

        end


        function [x] = solve(this, prdsid, A, b, varargin)
        % SOLVE solves a linear system with PARDISO
        %
        %    X = SOLVE(PRDS, PRDSID, A, B) solves the linear system A*X = B
        %    using the pardiso instance indicated by PRDSID. PRDS must be a
        %    pardiso object created calling the constructor of the class. The
        %    factorization of the matrix must be computed in advance calling
        %    the function FACTORIZE(PRDS, PRDSID, ...).
        % 
        %    X = SOLVE(PRDS, PRDSID, A, B, MNUM) solves the linear system 
        %    A*X = B for the matrix A relative to matrix number MNUM.
        %
        %    In case of symmetric matrices (matrix type = 2, -2, 6) only the
        %    lower triangular part of A should be passed as an argument.
        %
        %    See also PARDISO.PARDISO, PARDISO.FACTORIZE.
        %

            if nargin < 4
                error('Input must be the ParDiSo ID, one sparse matrix, and the RHS [and the matrix number].');
            end
            mnum = 1;
            if nargin == 5
                mnum = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y)  && size(A,1) == size(A,2));
            addRequired(p, 'b',      @(y) isnumeric(y) && size(y,1) == size(A,1));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, b, mnum);


            nrhs = size(b,2);
            perm = 0;
            err  = 0;

            if ~isreal(A)
                if isreal(b)
                    b = complex(b);
                end
            end

            %if issymmetric(A)
            %    A = tril(A);
            %end

            phase = 33;

            [~, this.iparm, ~, x, err, this.dparm] = ...
            callpardiso(prdsid, mnum, phase, ...
                        A, ...
                        perm, nrhs, this.iparm, this.msglvl, ...
                        b, 0, err, this.dparm);
            this.chkerr(err);


        end

        function [Ainv] = invert(this, prdsid, A, varargin)
        % INVERT performs selected inversion of a matrix with PARDISO
        %
        %    AINV = INVERT(PRDS, PRDSID, A) computes the selected inverse of the
        %    sparse matrix A on the ParDiSo instance indicated by PRDSID. PRDS
        %    must be a pardiso object created calling the constructor of the
        %    class. The factorization of the matrix must be computed in advance
        %    calling the function FACTORIZE(PRDS, PRDSID, ...).
        %    
        %    AINV = INVERT(PRDS, PRDSID, A, MNUM) computes the selected inverse
        %    for the matrix A relative to matrix number MNUM.
        %
        %    The selected inverse AINV is returned as a sparse matrix such that
        %    its sparsity pattern is the same as the one of the original matrix.
        %    In case of symmetric matrices (matrix type = 2, -2, 6) only the
        %    lower triangular part of A should be passed as an argument.
        % 
        %    See also PARDISO.PARDISO, PARDISO.FACTORIZE.
        %

            if nargin < 3
                error('Input must be the ParDiSo ID and one sparse matrix [and the matrix number].');
            end

            mnum = 1;
            if nargin == 4
                mnum = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y) && size(A,1) == size(A,2));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, mnum);

            err  = 0;
            perm = 0;
            dumm = [];
            nrhs = 1;

            %if issymmetric(A)
            %    A = tril(A);
            %end
            
            this.iparm(36) = 1;
            phase = -22;

            [Ainv, this.iparm, ~, ~, err, this.dparm] = ...
            callpardiso(prdsid, mnum, phase, ...
                        A, ...
                        perm, nrhs, this.iparm, this.msglvl, ...
                        dumm, dumm, err, this.dparm);
            this.chkerr(err);

        end



        function SC = schur(this, prdsid, A, s, varargin)
        % SCHUR returns the Schur-complement of A with respect to a subblock
        %
        %    SC = SCHUR(PRDS, PRDSID, A, S) computes the Schur complement of
        %    square sparse matrix A with respect to the block of size S at
        %    indices [end-s+1:end, end-s+1:end]. S must be a scalar greater than
        %    1 and less than the size of A. PRDS must be a pardiso object 
        %    created calling the constructor of the class.
        %    
        %    SC = SCHUR(PRDS, PRDSID, A, S) computes the Schur complement for
        %    matrix A relative to matrix number MNUM.
        % 
        %    See also PARDISO.PARDISO.
        %

            if nargin < 4
                error('Input must be the ParDiSo ID, one sparse matrix, and the size of the Schur-complement [and the matrix number].');
            end

            mnum = 1;
            if nargin == 5
                mnum = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'A',      @(y) issparse(y)  && size(A,1) == size(A,2));
            addRequired(p, 's',      @(y) isnumeric(y) && isscalar(y) && y > 1 && y < size(A,1));
            addRequired(p, 'mnum',   @(y) isnumeric(y) && isscalar(y) && y > 0 && y < this.maxfact);
            parse(p, prdsid, A, s, mnum);

            err  = 0;
            perm = 0;
            dumm = [];
            nrhs = 1;

            phase = 12;
            this.iparm(38) = s;     % Number of rows/cols for Schur complement.

            [SC, this.iparm, err, this.dparm] = ...
            callpardiso(prdsid, mnum, phase, ...
                        A, ...
                        perm, nrhs, this.iparm, this.msglvl, ...
                        dumm, dumm, err, this.dparm, s);

            this.iparm(38) = 0;
            this.chkerr(err);
        
        end

        function this = free_(this, prdsid, varargin)
        % FREE_ releases the internal memory in PARDISO
        %
        %    FREE_(PRDS, PRDSID) releases the memory relative to PRDSID. PRDS
        %    must be a pardiso object created calling the constructor of the
        %    class. 
        %
        %    See also PARDISO.PARDISO.
        %

            if nargin < 2
                error('Input must be the ParDiSo ID [and the memory ID]. If memory ID does not exist or is 0, all memory is realeased');
            end

            memid = 0;
            if nargin == 3
                memid = varargin{1};
            end

            p = inputParser;
            addRequired(p, 'prdsid', @(y) isnumeric(y) && isscalar(y) && y > 0);
            addRequired(p, 'memid',  @(y) isnumeric(y) && isscalar(y) && y >= 0 && y < this.maxfact);
            parse(p, prdsid, memid);

            err   = 0;
            dzero = 0;
            done_ = 1;
            dumm  = [];

            [err] = callpardiso(prdsid, memid);
            this.chkerr(err);
        end

        function ip = getiparm(this)
        % GETIPARM returns the IPARM array
        %
        %    IP = GETIPARM(PRDS) returns in array IP the column array PRDS.IPARM.
        %
        %    See also PARDISO.PARDISO.
        %

            ip = this.iparm;

        end

        function printiparm(this)
        % PRINTIPARM prints the IPARM array
        %
        %    PRINTIPARM(PRDS) prints out the column array PRDS.IPARM.
        %
        %    See also PARDISO.PARDISO.
        %

            fprintf('   iparm:\n');
            for k = 1:64
                fprintf('   %2d: %-10d\n', k, this.iparm(k));
            end
            fprintf('\n');
        end

        function dp = getdparm(this)
        % GETDPARM returns the DPARM array
        %
        %    DP = GETDPARM(PRDS) returns in array DP the column array PRDS.DPARM.
        %
        %    See also PARDISO.PARDISO.
        %

            dp = this.dparm;

        end

        function printdparm(this)
        % PRINTDPARM prints the IPARM array
        %
        %    PRINTDPARM(PRDS) prints out the column array PRDS.DPARM.
        %
        %    See also PARDISO.PARDISO.
        %

            fprintf('   dparm:\n');
            for k = 1:64
                fprintf('   %2d: %-0.4ed\n', k, this.dparm(k));
            end
            fprintf('\n');
        end


        
        function pc = getcounter(this)
        % GETCOUNTER returns the number of ParDiSo instances currently available
        %
        %    PC = GETCOUNTER(PRDS) returns the counter indicating the number of
        %    ParDiSo instances currently allocated in memory. PC indicates the
        %    instances allocated but do not take trace of the ones released by
        %    the method PARDISO.FREE_().
        %
        %    See also PARDISO.PARDISO, PARDISO.FREE_.
        %

            pc = this.prds_count;
        end
        

    end   % end of public methods


    methods (Static)

        function s = isstructsymm(A)
        % ISSTRUCTSYMM checks if a matrix is structurally symmetric
        %
        %    S = ISSTRUCTSYMM(A) returns true or false whether A is structurally
        %    symmetric or not. A must be a sparse, square matrix.
        %
        %    See also PARDISO.PARDISO, PARDISO.MAKESTRUCTSYMM.
        %

            if ~issparse(A) || size(A,1) ~= size(A,2)
                error('isstructsymm:invalidInputType', ...
                    'Input must be a square, sparse matrix.');
            end

            s = issymmetric(spones(A));
    
        end

        function S = makestructsymm(A)
        % MAKESTRUCTSYMM returns a structurally symmetric version of the input
        %
        %    S = MAKESTRUCTSYMM(A) returns matrix A added with the minimum
        %    number of nonzeros such that it becomes structurally symmetric. The
        %    nonzero values added to obtain a symmetric sparsity pattern are all
        %    equal to REALMIN.
        %
        %    See also PARDISO.PARDISO, PARDISO.ISSTRUCTSYMM, REALMIN.
        %

            D = spones(A) - spones(A.');
            D = (D == -1);
            S = A + realmin * D;
            fprintf('Warning: Adding %d nonzeros (%0.1e%%) to make the matrix structurally symmetric.\n', nnz(D), nnz(D)/nnz(A)*100);

        end

	function mtype = getpardisomtype(A)
        % GETPARDISOMTYPE returns the ParDiSo matrix type
        %
        %    MTYPE = GETPARDISOMTYPE(A) returns the ParDiSo matrix type of the
        %    input matrix A.
        %
        %    See also PARDISO.PARDISO.
        %

	    if ~issparse(A) || size(A,1) ~= size(A,2)
		error('isstructsymm:invalidInputType', ...
		    'Input must be a square, sparse matrix.');
	    end

	    if isreal(A)
                if issymmetric(A)
                    mtype = -2;
                elseif pardiso.isstructsymm(A)
                    mtype = 1;
                else
                    mtype = 1;  % Matrix will be modified to be struct-symm
                end
	    else
                if issymmetric(A)
                    mtype = 6;
                elseif pardiso.isstructsymm(A)
                    mtype = 3;
                elseif ishermitian(A)
                    mtype = -4;
                else
                    mtype = 3;  % Matrix will be modified to be struct-symm
                end
	    end


	end


    end   % end of static methods

    methods (Access = protected)

        function chkerr(this, error_)

            if error_ ~= 0

                fprintf('ParDiSo error: %2d --- ', error_);

                if error_ == -1
                    fprintf('Input inconsistent\n');
                elseif error_ == -2
                    fprintf('Not enough memory\n');
                elseif error_ == -3
                    fprintf('Reordering problem\n');
                elseif error_ == -4
                    fprintf('Zero pivot, numerical factorization or iterative refinement problem\n');
                elseif error_ == -5
                    fprintf('Unclassified (internal) error\n');
                elseif error_ == -6
                    fprintf('Preordering failed (matrix type 11, 13 only)\n');
                elseif error_ == -7
                    fprintf('Diagonal matrix problem\n');
                elseif error_ == -8
                    fprintf('32-bit integer overflow problem\n');
                elseif error_ == -10
                    fprintf('No license file pardiso.lic found\n');
                elseif error_ == -11
                    fprintf('License is expired.\n');
                elseif error_ == -12
                    fprintf('Wrong username or hostname.\n');
                end

                exit;

            else
                    % Nothing!
            end

        end
    end  % end of protected methods

end




