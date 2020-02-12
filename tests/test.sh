# Use "&&" so the test will fail if either python call fails. 
cd tests/system && python -m unittest && cd .. && python -m unittest && cd ..

