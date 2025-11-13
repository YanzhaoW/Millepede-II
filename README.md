# Millepede-II

Millepede II is a package for linear least squares fits with a large number of parameters. Developed for the alignment and calibration of tracking detectors.

## CMake configuration

Use the git submodule to include this repository into your project:

```bash
git submodule add https://github.com/YanzhaoW/Millepede-II.git millepede
```

```cmake
add_subdirectory(millepede)

target_link_libraries(main PRIVATE millepede::mille)
```

## Usage of Mille interface

The main classes of the Mille interface:

- `MilleWriter`: for writing the data to the binary file, as the input for the pede program
- `ResultReader`: for reading the output parameter from pede program
- `SteerWriter`: for writing the configuration file for the pede program

### MilleWriter

`MilleWriter` class is the main interface class to write the data to a binary file.

Its public interfaces can be seen below:

```cpp
namespace millepede
{
    class MilleWriter
    {
        explicit MilleWriter(const std::string& out_filename, bool is_binary = true, bool is_write_zero = false);
        void set_buffer_size(std::size_t buffer_size);
        void mille(const MilleDataPoint& data_point);
        void special(const std::vector<std::pair<int, float>>& special_data);
        void end();
        void close();
    };
}
```

#### Constructor

```cpp
explicit MilleWriter(const std::string& out_filename, bool is_binary = true, bool is_write_zero = false);
```

The constructor sets the filename of the output binary file. If the `is_binary` parameter is true, the output data will be saved as the binary format. If the `is_write_zero` parameter is false, only non-zero values from the member function call `mille` are written to the output file.

### mille

```cpp
void mille(const MilleDataPoint& data_point);
```

The `mille` member function adds the data point `data_point` to the current entry. Data point has the type `MilleDataPoint` and has the following member variables:

```cpp
struct MilleDataPoint
{
    std::vector<float> locals;                  // local derivatives
    std::vector<std::pair<int, float>> globals; // global label and derivatives pair
    float measurement = 0.;                     // measurement corresponding to the error value
    float sigma = 1.;                           // error value
};
```

- `locals`: a vector of first-order derivatives of local variables
- `globals`: a vector of global variable related pairs, with the first element to be the global parameter ID and second element to be the corresponding first-order derivative
- `measurement`: measurement of the data point, which is just a constant term without local and global parameters.
- `sigma`: the error corresponding to the measurement.

### end

```cpp
void end();
```

This saves the current entry to the output file. It also resets all internal buffers to be ready for the next entry.

### close

```cpp
void close();
```

This manually calls the `close` member function of the internal file ostream of type `std::ofstream`. This member function will be automatically called when the `MilleWriter` object is deleted.

### ResultReader

This class can be used to read the result file from the pede program, which is name `millepede.res` by default.

Public interfaces:

```cpp
    class ResultReader
    {
      public:
        ResultReader() = default;
        void set_filename(std::string_view filename);

        void read();
        void print();
        [[nodiscard]] auto get_pars() const -> const std::unordered_map<int, ParResultEntry>&;
    };
```

#### Constructor

```cpp
ResultReader() = default;
```

Default constructor set the filename to be an empty string.

#### set_filename

```cpp
void set_filename(std::string_view filename);
```

This sets the filename of the resulting parameter output from the pede program.

#### read

```cpp
void read();
```

This function triggers the reading the parameter from the file.

#### get_pars

```cpp
[[nodiscard]] auto get_pars() const -> const std::unordered_map<int, ParResultEntry>& { return par_results_; }
```

This function returns a constant reference to the parameters that has been read from the file. The parameters are stored in a hash table with the type `std::unordered_map<int, ParResultEntry>`. The type of the value is:

```cpp
    struct ParResultEntry
    {
        int par_id = 0;
        float value = 0.F;
        float pre_sigma = 0.F;
        float value_diff = 0.F;
        float error = 0.F;
    };
```

where

- `par_id`: the ID the global parameter.
- `value`: the value the global parameter.
- `pre_sigma`: the user provided sigma value for the global parameter.
- `value_diff`: the value difference between the user provided value and the fitted value.
- `error`: the fitted error value of the global parameter.

### SteerWriter

`SteerWriter` class can be used to create the configuration file for the pede program

Public interfaces:

```cpp
    class SteerWriter
    {
      public:
        SteerWriter() = default;

        // setters:
        void set_filepath(std::string_view filepath);
        void set_data_filepath(std::string_view filepath);
        void set_parameter_file(std::string_view filename);
        void set_working_dir(std::string_view dir);

        void add_parameter_default(int par_id, const std::pair<float, float>& values);
        void add_method(PedeMethod method, const std::pair<float, float>& values);
        void add_other_options(std::vector<std::string> options);

        void write();
    };
```

#### Constructor

```cpp
SteerWriter() = default;
```

The default constructor set the default values:

- `filepath`: "steer.txt"
- `data_filepath`: "data.bin"
- `parameter_file`: "parameters.txt"
- `working_dir`: "." (current folder)

#### add_parameter_default

```cpp
void add_parameter_default(int par_id, const std::pair<float, float>& values);
```

This adds the default values for the global parameters with the parameter ID and a pair of values, with the first one to be the default value and the second to the error value.

#### add_method

```cpp
void add_method(PedeMethod method, const std::pair<float, float>& values);
```

This sets the fitting algorithm for the pede program. The method has the enum type `PedeMethod`, which is defined as:

```cpp
enum class PedeMethod : uint8_t
{
    inversion,
    diagonalization,
    fullGORES,
    sparseGMRES,
    cholesky,
    bandchooseby,
    HIP
};
```

The values from each method can be looked up in the millepede-II manual page.

#### add_other_options

```cpp
void add_other_options(std::vector<std::string> options);
```

Add other options with a vector of strings. These strings are written into the configuration file in a single row separated by the space.

#### write

```cpp
void write();
```

This writes the configurations to the file. It should be called only after all the parameters and the method are set.

