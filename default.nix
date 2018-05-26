{ python3Packages }:

python3Packages.buildPythonPackage {
  pname = "photon-tools";
  version = "0.1";
  src = ./.;
  propagatedBuildInputs = with python3Packages; [numpy scipy matplotlib squmfit];
  nativeBuildInputs = with python3Packages; [cython];
}
