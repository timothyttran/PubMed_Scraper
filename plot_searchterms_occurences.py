"""
Author: Timothy Tran
Date: 9/18/2024

This script creates:
- a bar chart showcases the top/bottom search terms by average growth rate over a selected time period
- top_terms_growth_rates.csv : A list of the search terms with highest average growth rate

Note that due to the large number of search terms (>100), we restrict only a certain amount to show in the diagram.
Hence the "selected_terms" global variable (can use "avoided_terms" if we want to do all terms except for some).
"""
import pandas as pd
import matplotlib.pyplot as plt

top_and_bottom_N_terms = 15 # Show the top N and bottom N term in the bar graph

# Read the CSV file
df = pd.read_csv("occurrences_by_popularity.csv").iloc[1:, 1:]

selected_terms = {"Social anxiety", "Major depressive disorder", "Neuroimaging",
                  "Clinical psychology", "Substance dependence", "School psychology",
                  "Alzheimer's disease"}

avoided_terms = {"Behavior", "Psychologist", "Social psychology", "Personality psychology", "Health psychology", "Cognition"}
avoided_terms = {"Pervasive developmental disorder not otherwise specified"}

# Initialize a dictionary to store data for each term
term_data = {}

# Iterate over each column (year) in the DataFrame
for year in df.columns[1:-2]:
    print(year)
    # Iterate over each row in the DataFrame for the corresponding year
    for term_occurence in df[year]:
        # Split the cell value into term and occurrence
        term, occurrence = term_occurence.strip().split(',')
        term = term.strip('"')
        occurrence = int(occurrence.strip())

        # Add the data to the dictionary
        if term not in avoided_terms:
            if term not in term_data:
                term_data[term] = {'years': [], 'occurrences': []}
            term_data[term]['years'].append(int(year))  # Convert year to int
            term_data[term]['occurrences'].append(occurrence)

# Calculate the average annual growth rate for each term
average_growth_rate = {}
for term, data in term_data.items():
    if len(data['occurrences']) > 1:  # Ensure there are at least two data points
        values = data['occurrences']
        growth_rates = [(values[i] - values[i-1]) / values[i-1] for i in range(1, len(values))]
        average_growth_rate[term] = sum(growth_rates) / len(growth_rates) * 100

# Select the top N terms based on the average growth rate
top_terms = sorted(average_growth_rate, key=average_growth_rate.get, reverse=True)
top_and_bottom_terms = top_terms[:top_and_bottom_N_terms] + top_terms[-top_and_bottom_N_terms:]
rates = [average_growth_rate[term] for term in top_and_bottom_terms]

# Create CSV file for a table of the top search terms and their growth rates
rates_top_terms = [round(average_growth_rate[term], 2) for term in top_terms]
df = pd.DataFrame({
    'Term': top_terms,
    'Growth Rate (%)': rates_top_terms
})
# Add a Rank column
df['Rank'] = range(1, len(df) + 1)
# Reorder columns to have Rank first
df = df[['Rank', 'Term', 'Growth Rate (%)']]
# Write the DataFrame to a CSV file
df.to_csv('top_terms_growth_rates.csv', index=False)

# Plot the data for each term
for term, data in term_data.items():
    plt.plot(data['years'], data['occurrences'], label=term)

# Add labels and legend
plt.xlabel('Year')
plt.ylabel('Occurrence')
plt.title('Search Term Occurrence Over Years')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()

# Plot the growth rate for the top terms in the second window
plt.figure(figsize=(12, 8))
plt.bar(top_and_bottom_terms, rates, color='blue')
plt.ylabel('Growth Rate (%)')
plt.xlabel('Search Term / Article Name')
plt.title(f'Top and Bottom {top_and_bottom_N_terms} Search Terms by Average Annual Growth Rate')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.grid(True, axis='y')
plt.show()