from mosden.utils.csv_handler import CSVHandler
from mosden.utils.input_handler import InputHandler
import os
import pytest

def test_csv_handler():
    """
    Test the CSVHandler class for reading and writing CSV files.
    """
    test_path = './tests/unit/output/test.csv'
    csv_handler = CSVHandler(test_path, overwrite=True)

    # Test writing to a CSV file
    data = {'Xe135': {'value': 1.0, 'uncertainty': 0.1},
            'U235': {'value': 100.0, 'uncertainty': 5.0}}
    csv_handler.write_csv(data)
    
    # Test reading from the same CSV file
    read_data = csv_handler.read_csv()
    
    assert read_data == data, f"Expected {data}, but got {read_data}"
 
    return

def test_input_handler():
    """
    Test the InputHandler class for reading input files.
    """
    input_path = 'test_input.json'
    
    # Create a test input file
    with open(input_path, 'w') as f:
        f.write('{"key": "value"}')
    
    input_handler = InputHandler(input_path)
    
    # Test reading the input file
    data = input_handler.read_input(check=False)
    
    assert data == {"key": "value"}, f"Expected {{'key': 'value'}}, but got {data}"

    # Test exceptions in input handling
    with pytest.raises(KeyError):
        data = input_handler.read_input(check=True)

    os.remove(input_path)
    
    return