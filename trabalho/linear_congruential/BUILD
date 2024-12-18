cc_binary(
    name = "main",
    srcs = ["main.cc"],
    copts = [
        "-I/usr/include/eigen3",
        # "-I/usr/include",  # Update this path to where Matplot++ headers are located
        "-std=c++17",
        "-g",
    ],
    linkopts = [
        # "-L/usr/lib",  # Update this path to where Matplot++ library is located
        "-lmatplot", # Update this to the actual Matplot++ library name if different

    ],
)