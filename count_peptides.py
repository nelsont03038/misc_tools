grades = ['b', 'b', 'f', 'c', 'b', 'a', 'a', 'd', 'c', 'd', 'a', 'a', 'b']

print(list(map(lambda letter: {letter: grades.count(letter)}, set(grades))))

for item in occurrence:
    print(f"{item[0]}={item[1]}")
	





from collections import defaultdict

grades = ['b', 'b', 'f', 'c', 'b', 'a', 'a', 'd', 'c', 'd', 'a', 'a', 'b']
occurrence = defaultdict(lambda: 0)

for character in grades:
    occurrence[character] += 1

for key, value in occurrence.items():
    print(f"{key}={value}")




high_school_semester_grades = ["B+", "A", "B+", "A", "A+", "A-", "A"]
from collections import Counter
from collections import defaultdict


# needs no imports
dict_count_occurences_of_letters = {}

for letter in high_school_semester_grades:
    if letter in dict_count_occurences_of_letters:
        dict_count_occurences_of_letters[letter] += 1
    else:
        dict_count_occurences_of_letters[letter] = 1
sorted(dict_count_occurences_of_letters.items(), key=lambda x: x[1], reverse=True)



defaultdict_count_occurences_of_letters = defaultdict(int)
for letter in high_school_semester_grades:
    defaultdict_count_occurences_of_letters[letter] += 1
sorted(defaultdict_count_occurences_of_letters.items(), key=lambda x: x[1], reverse=True)

	


counter_count_of_letters = Counter(high_school_semester_grades)
counter_count_of_letters.most_common()


