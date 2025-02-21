import ipywidgets as widgets

# Global stored values dictionary
stored_values = {
    'Simulation type': None,
    'Sequence': None,
    'Sequence 1': None,
    'Sequence 2': None,
    'Ratio': None,
    'Processed_Ratio': None,
    'peptide_pairs': None
}

# Create widgets
checkbox_a = widgets.Checkbox(
    value=False,
    description='Single peptide',
    disabled=False
)

checkbox_b = widgets.Checkbox(
    value=False,
    description='Co-assembly',
    disabled=False
)

text_input_a = widgets.Text(
    description='Sequence:',
    disabled=False,
    layout={'visibility': 'hidden', 'width': '800px'},
    placeholder='Enter sequence (up to 100 characters)'
)

text_input_b1 = widgets.Text(
    description='Sequence 1:',
    disabled=False,
    layout={'visibility': 'hidden', 'width': '800px'},
    placeholder='Enter sequence 1 (up to 100 characters)'
)

text_input_b2 = widgets.Text(
    description='Sequence 2:',
    disabled=False,
    layout={'visibility': 'hidden', 'width': '800px'},
    placeholder='Enter sequence 2 (up to 100 characters)'
)

text_input_b3 = widgets.Text(
    description='Ratio:',
    disabled=False,
    layout={'visibility': 'hidden', 'width': '800px'},
    placeholder='Enter ratio'
)

def process_ratio(ratio_str):
    try:
        numbers = ratio_str.replace(' ', '').split('-')
        if len(numbers) == 2:
            num1 = int(numbers[0]) * 100
            num2 = int(numbers[1]) * 100
            return f"{num1}-{num2}"
    except:
        return ratio_str
    return ratio_str

def update_stored_values():
    if checkbox_a.value:
        stored_values['Simulation type'] = 'Single peptide'
        stored_values['Sequence'] = text_input_a.value
        stored_values['Sequence 1'] = text_input_a.value
        stored_values['Sequence 2'] = 'A'
        stored_values['Ratio'] = '10-0'
        stored_values['Processed_Ratio'] = process_ratio('10-0')
        if text_input_a.value:
            stored_values['peptide_pairs'] = f"{text_input_a.value}-A"
    elif checkbox_b.value:
        stored_values['Simulation type'] = 'Co-assembly'
        stored_values['Sequence'] = None
        stored_values['Sequence 1'] = text_input_b1.value
        stored_values['Sequence 2'] = text_input_b2.value
        stored_values['Ratio'] = text_input_b3.value
        stored_values['Processed_Ratio'] = process_ratio(text_input_b3.value)
        if text_input_b1.value and text_input_b2.value:
            stored_values['peptide_pairs'] = f"{text_input_b1.value}-{text_input_b2.value}"

def on_checkbox_change(change):
    if change['owner'] == checkbox_a and change['new']:
        checkbox_b.value = False
        text_input_a.layout.visibility = 'visible'
        text_input_b1.layout.visibility = 'hidden'
        text_input_b2.layout.visibility = 'hidden'
        text_input_b3.layout.visibility = 'hidden'
    elif change['owner'] == checkbox_b and change['new']:
        checkbox_a.value = False
        text_input_a.layout.visibility = 'hidden'
        text_input_b1.layout.visibility = 'visible'
        text_input_b2.layout.visibility = 'visible'
        text_input_b3.layout.visibility = 'visible'
    update_stored_values()

def on_text_change(change):
    update_stored_values()

def create_widget_layout():
    widget = widgets.VBox([
        checkbox_a,
        checkbox_b,
        text_input_a,
        text_input_b1,
        text_input_b2,
        text_input_b3
    ])
    return widget

# Register callbacks
def initialize_callbacks():
    checkbox_a.observe(on_checkbox_change, names='value')
    checkbox_b.observe(on_checkbox_change, names='value')
    text_input_a.observe(on_text_change, names='value')
    text_input_b1.observe(on_text_change, names='value')
    text_input_b2.observe(on_text_change, names='value')
    text_input_b3.observe(on_text_change, names='value')