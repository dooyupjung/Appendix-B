#define NAME 200
#define SPEC 2000
#define PROT 2000

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <direct.h>

int main(void)
{
	FILE *fin, *fout;
	char seq_name[NAME], new_name[NAME], **species, newdir[NAME];
	int i, j, start = 0, end = 0, len, num_species = 0, flag, num_proteins[SPEC] = { 0, }, total_proteins = 0;
	int **protein_index, cnt, k, l, m;

	species = (char **)malloc(sizeof(char*) * SPEC);
	protein_index = (int **)malloc(sizeof(int*) * SPEC);
	for (i = 0; i < SPEC; i++) {
		species[i] = (char *)malloc(sizeof(char) * NAME);
		protein_index[i] = (int *)malloc(sizeof(int) * PROT);
	}

	// step 1: species�� ������ ���� �ľ�
	fopen_s(&fin, "vertebrates.fasta", "r");
	printf("Open!\n");

	while (!feof(fin)) {
		//�ʱ�ȭ �� �Է� �ޱ�
		for (i = 0; i < NAME; i++)
			seq_name[i] = '\0';

		fgets(seq_name, NAME, fin);

		if (seq_name[0] == '>') {
			total_proteins++;

			// ���� �ʱ�ȭ
			for (i = 0; i < NAME; i++)
				new_name[i] = '\0';
			
			// species �̸� ����
			len = strlen(seq_name);
			for (i = 1; i < len; i++) {
				if (seq_name[i] == '[')	start = i;
				if (seq_name[i] == ']') end = i;
			}

			for (i = 0; i < end - start - 1; i++)
				new_name[i] = seq_name[i + start + 1];

			new_name[i] = '.';
			new_name[i + 1] = 'f';
			new_name[i + 2] = 'a';
			new_name[i + 3] = 's';
			new_name[i + 4] = 't';
			new_name[i + 5] = 'a';

			// Species �̸� �������� �� ����ó��!
			for (i = 0; i < NAME; i++)
				if (new_name[i] == ':' || new_name[i] == '/')
					new_name[i] = ' ';

			// new species���� Ȯ�� : flag�� 1�̸� new, 0�̸� �̹� ����
			flag = 1;

			for (i = 0; i < num_species; i++) {
				if (strcmp(new_name, species[i]) == 0) {
					flag = 0;
					protein_index[i][num_proteins[i]] = total_proteins;
					num_proteins[i]++;
					break;
				}
			}

			// new species�� ���, ��Ͽ� add-up
			if(flag == 1) {
				strcpy_s(species[num_species], NAME, new_name);
				protein_index[num_species][0] = total_proteins;
				num_proteins[num_species]++;
				num_species++;
			}
		}
	}

	fclose(fin);

	system("mkdir Fasta_species\\");
	fopen_s(&fout, "Fasta_species/summary.txt", "w");
	
	fprintf(fout, "Split summary of FASTA file:\n");
	fprintf(fout, "Total %d proteins were splitted according to %d species!\n\n", total_proteins, num_species);
	for (i = 0; i < num_species; i++) {
		fprintf(fout, "* %s: %d proteins\n", species[i], num_proteins[i]);
		fprintf(fout, "\t>Index number:");
		for (j = 0; j < num_proteins[i]; j++)
			fprintf(fout, " %d", protein_index[i][j]);
		fprintf(fout, "\n");
	}
	
	fclose(fout);

	fopen_s(&fout, "Fasta_species/filename.txt", "w");
	fprintf(fout, "%s\n%d", species[0], num_proteins[0]);
	for (i = 1; i < num_species; i++)
		fprintf(fout, "\n%s\n%d", species[i], num_proteins[i]);
	fclose(fout);

	// step 2: FASTA file splitting according to species
	
	for (i = 0; i < num_species; i++) {
		for (m = 0; m < NAME; m++)
			newdir[m] = '\0';

		strcat_s(newdir, NAME, "Fasta_species/");
		strcat_s(newdir, NAME, species[i]);
		fopen_s(&fout, newdir, "w");

		fopen_s(&fin, "vertebrates.fasta", "r");
		flag = 0;
		cnt = 0;
		l = 0;
		while (!feof(fin)) {
			//�ʱ�ȭ �� �Է� �ޱ�
			for (j = 0; j < NAME; j++)
				seq_name[j] = '\0';

			fgets(seq_name, NAME, fin);

			if (seq_name[0] == 10)	// blank line�� ����
				continue;

			else {
				for (k = 1; k < NAME; k++) {
					if (seq_name[k] == 10) {
						seq_name[k] = '\0';
						break;
					}
				}

				if (seq_name[0] == '>') {	// start line�� ���
					cnt++;
					if (cnt == protein_index[i][l]) {		// index�� ��ġ�ϸ�,
						if (l >= 1)
							fprintf(fout, "\n");

						fprintf(fout, "%s\n", seq_name);
						l++;
						flag = 1;			// general line ��� on!
					}

					else {  // index�� ��ġ���� ������,
						flag = 0;	// general line ��� off!
					}

				}

				else {						// general line�� ���
					if (flag == 1) {		// ��� on ���¶��
						fprintf(fout, "%s\n", seq_name);	// ���!
					}
				}
			}
		}

		fclose(fin);
		fclose(fout);
	}

	for (i = 0; i < SPEC; i++) {
		free(species[i]);
		free(protein_index[i]);
	}
	free(species);
	free(protein_index);
}