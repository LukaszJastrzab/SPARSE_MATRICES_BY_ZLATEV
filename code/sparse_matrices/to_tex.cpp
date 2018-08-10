//template <typename TYPE>
//class input_storage_scheme
//{
//public:
//    /// number of added elements - sizes of arrays (std::vector): AORIG, RNORIG and CNORIG
//    const size_t NNZ;
//    /// assumed sizes of stored matrix
//    const size_t number_of_rows, number_of_columns;
//
//private:
//    /// array of ORIGINAL values of the input matrix A
//    std::vector<TYPE> AORIG;
//    /// array of ORIGINAL row numbers of the input matrix A (indexed from 0)
//    std::vector<int> RNORIG;
//    /// array of ORIGINAL column numbers of the input matrix A (indexed from 0)
//    std::vector<int> CNORIG;
//
//    /// disable default construcor
//    input_storage_scheme();
//
//public:
//    /// main constructor - requires only sizes of matrix
//    input_storage_scheme(int _number_of_rows, int _number_of_columns);
//    /// method to adding elements
//    void add_element(TYPE value, size_t row, size_t col);
//
//    /// other public methods...
//private:
//    /// private methods...
//
//};

//template <typename TYPE>
//class dynamic_storage_scheme
//{
//private:
//    //======== MATRIX - BASIC INFORMATIONS ========
//    const size_t number_of_rows, number_of_columns;
//
//    //========================= Row-Ordered List (ROL) =========================
//    size_t NROL;    /// size of row-ordered list
//    size_t LROL;    /// last not free position in ROL
//    size_t CROL;    /// number of non-zeros actualy stored in ROL
//    TYPE *ALU;      /// array of values of the elements stored in scheme
//    int *CNLU;      /// array of column numbers of the elements stored in scheme
//
//    //====================== Column-Ordered List (COL) ======================
//    size_t NCOL;    /// size of column-ordered list
//    size_t LCOL;    /// last not free position in COL
//    size_t CCOL;    /// number of non-zeros actualy stored in COL
//    int *RNLU;      /// array of row numbers of the elements stored in scheme
//
//    //================= INTEGRITY ARRAYS =================
//    size_t NHA;     /// size of integrity arrays
//    TYPE *PIVOT;    /// table of pivots
//    int **HA;       /// Array of pointers and permutations
//
//    std::string logged_errors;      /// string used to logging errors
//    dynamic_storage_scheme();   /// disable default constructor
//
//public:
//    /// Constructor - input_storage_scheme is required and two floats
//    ///               that determines sizes of storage lists in flollowing way:
//    ///               NROL = mult1 * ISS->NNZ, NCOL = mult2 * NROL
//    dynamic_storage_scheme(const input_storage_scheme<TYPE>& ISS,
//                           double mult1, double mult2 = 0.7);
//    /// Destructor
//    ~dynamic_storage_scheme();
//
//    /// other public methods...
//private:
//    /// private methods...
//};

//template <typename TYPE>
//class dynamic_storage_scheme
//{
//private:
//    //======== MATRIX - BASIC INFORMATIONS ========
//    /// sizes of stored matrix
//    const size_t number_of_rows, number_of_columns;
//
//    //============== ROW-ORDERED LIST (ROL) =============
//    size_t NROL;    /// size of row-ordered list
//    size_t LROL;    /// last not free position in ROL
//    size_t CROL;    /// number of non-zeros actualy stored in ROL
//    TYPE *ALU;      /// array of values of the elements stored in scheme
//    int *CNLU;      /// array of column numbers of the elements stored in scheme
//
//    //=============== COLUMN-ORDERED LIST (COL) ===============
//    size_t NCOL;    /// size of column-ordered list
//    size_t LCOL;    /// last not free position in COL
//    size_t CCOL;    /// number of non-zeros actualy stored in COL
//    int *RNLU;      /// array of row numbers of the elements stored in scheme
//
//    //============== INTEGRITY ARRAYS =============
//    size_t NHA;     /// size of integrity arrays
//    TYPE *PIVOT;    /// place for important values
//    int **HA;       /// array of pointers and permutations
//
//    //============== OTHER FIELDS ==============
//    DYNAMIC_STATE dynamic_state;   /// flag indicating current state of the scheme
//    double omega;                  /// relaxation parameter for SOR method
//    std::string logged_errors;     /// output string used to logging errors
//
//    /// Disable default constructor
//    dynamic_storage_scheme();
//
//public:
//    /// Constructor - input_storage_scheme and two floats are required to determine 
//    /// lists sizes in flollowing way: NROL = mult1 * ISS->NNZ, NCOL = mult2 * NROL
//    dynamic_storage_scheme(const input_storage_scheme<TYPE>& ISS,
//                           double mult1, double mult2 = 0.7);
//    /// Destructor
//    ~dynamic_storage_scheme();
//    // other public methods...
//private:
//    // private methods...
//};

//template <typename TYPE>
//class input_storage_scheme
//{
//public:
//    /// number of added elements - sizes of arrays (std::vector): AORIG, RNORIG and CNORIG
//    const size_t NNZ;
//    /// assumed by bellow variables order of a matrix
//    const size_t order;
//    /// assumed sizes of stored matrix
//    const size_t number_of_rows, number_of_columns;
//
//private:
//    /// array of ORIGinal values of the input matrix A
//    std::vector<TYPE> AORIG;
//    /// array of ORIGinal row numbers of the input matrix A (indexed from 0)
//    std::vector<int> RNORIG;
//    /// array of ORIGinal column numbers of the input matrix A (indexed from 0)
//    std::vector<int> CNORIG;
//
//public:
//    /// default constructor
//    input_storage_scheme(TYPE val)
//        :
//        NNZ(1),
//        order(1),
//        number_of_rows(1),
//        number_of_columns(1)
//    {
//        AORIG.push_back(val);
//        RNORIG.push_back(0);
//        CNORIG.push_back(0);
//    }
//    /// constructor which requires sizes of matrix
//    input_storage_scheme(int _number_of_rows, int _number_of_columns)
//        :
//        number_of_rows(_number_of_rows),
//        number_of_columns(_number_of_columns),
//        order(_number_of_rows < _number_of_columns ? _number_of_rows : _number_of_columns),
//        NNZ(0)
//    {
//    }
//    /// method to adding elements
//    void add_element(TYPE value, size_t row, size_t col);
//
//    /// other public methods...
//private:
//    /// private methods...
//};