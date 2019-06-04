#define _CRT_SECURE_NO_WARNINGS

#define NAME 200
#define MAX 15000
#define LENGTH 3
#define AATYPE 26

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <direct.h>

int main(void) {
	FILE *fin, *fseq, *fout, *ffile1, *ffile2, *fsum1, *fsum2, *fspace;
	char filename[NAME], seq[NAME], newdir[NAME], whole[MAX], cmp1[LENGTH + 1], cmp2[LENGTH + 1];
	int i, j, k, len, m, n, end = 0, p, q, r, motif_type, flag, num_proteins, temp;
	int motif_db[LENGTH][AATYPE][AATYPE];

	cmp1[LENGTH] = '\0';
	cmp2[LENGTH] = '\0';

	// Species ��� �����ִ� ���� ����
	fin = fopen("Fasta_species/filename.txt", "r");

	// Raw_motifs ���� �����
	system("mkdir Raw_motifs\\");

	// Raw_motifs/filename.txt, summary.txt �����
	ffile1 = fopen("Raw_motifs/filename.txt", "w");
	fsum1 = fopen("Raw_motifs/summary.txt", "w");

	// Without_space ���� �����
	system("mkdir Without_space\\");

	// Without_space/filename.txt �����
	ffile2 = fopen("Without_space/filename.txt", "w");
	fsum2 = fopen("Without_space/summary.txt", "w");

	fprintf(fsum1, "Summary of motif counting:\n");
	fprintf(fsum2, "Summary of sequence extraction without space:\n");

	while (!feof(fin)) {
		fprintf(fsum1, "\n");
		fprintf(fsum2, "\n");

		for (i = 0; i < NAME; i++) {
			filename[i] = '\0';
			newdir[i] = '\0';
		}

		fgets(filename, NAME, fin);
		fscanf(fin, "%d\n", &num_proteins);

		// ���๮�� NULL�� ����ó��
		for (i = 0; i < NAME; i++) {
			if (filename[i] == 10) {
				filename[i] = '\0';
				break;
			}
		}

		// Species���� �������� FASTA���� �ϳ��� �ҷ�����
		strcat_s(newdir, NAME, "Fasta_species/");
		strcat_s(newdir, NAME, filename);
		fopen_s(&fseq, newdir, "r");

		// MOTIF ������ ��ϵ� (species ����).motif ���� ���� 
		filename[strlen(filename) - 5] = 'm';
		filename[strlen(filename) - 4] = 'o';
		filename[strlen(filename) - 3] = 't';
		filename[strlen(filename) - 2] = 'i';
		filename[strlen(filename) - 1] = 'f';

		for (i = 0; i < NAME; i++)
			newdir[i] = '\0';

		strcat_s(newdir, NAME, "Raw_motifs/");
		strcat_s(newdir, NAME, filename);

		fout = fopen(newdir, "w");

		// MOTIF ������ ��ϵ� (species ����).motif �����̸� ��� ����
		fprintf(ffile1, "%s\n%d\n", filename, num_proteins);

		// ������� ���� ������ ��ϵ� (species ����).space ���� ���� 
		filename[strlen(filename) - 5] = 's';
		filename[strlen(filename) - 4] = 'p';
		filename[strlen(filename) - 3] = 'a';
		filename[strlen(filename) - 2] = 'c';
		filename[strlen(filename) - 1] = 'e';

		for (i = 0; i < NAME; i++)
			newdir[i] = '\0';

		strcat_s(newdir, NAME, "Without_space/");
		strcat_s(newdir, NAME, filename);

		fspace = fopen(newdir, "w");

		// MOTIF ������ ��ϵ� (species ����).motif �����̸� ��� ����
		fprintf(ffile2, "%s\n%d\n", filename, num_proteins);
		
		flag = 0;	// ���๮�� ���� ���� ���

		while (!feof(fseq)) {		// Species���� �������� FASTA ���� ���� ����
			// �迭 �ʱ�ȭ �� ���ڿ� �� �� �Է¹ޱ�
			for (j = 0; j < NAME; j++)
				seq[j] = '\0';
			fgets(seq, NAME, fseq);
			
			// �Է¹��� ���ڿ� ���๮�� ����
			for (j = 0; j < NAME; j++) {
				if (seq[j] == 10) {
					seq[j] = '\0';
					break;
				}
			}

			if (seq[0] == '>') {		// �Է¹��� ���ڿ��� start line�� ���, �״�� MOTIF, SPACE ���Ͽ� ����ϰ� ���� �ʱ�ȭ
				printf("Treating %s...\n", seq);
				
				if (flag == 0)			// Seq ���̿� ���๮�� ����
					flag = 1;

				else {
					fprintf(fout, "\n");
					fprintf(fspace, "\n");
				}

				fprintf(fout, "%s\n", seq);
				fprintf(fspace, "%s\n", seq);

				fprintf(fsum1, "%s\n", &seq[1]);
				fprintf(fsum2, "%s\n", &seq[1]);

				len = 0;
				end = 0;
				for (k = 0; k < MAX; k++)
					whole[k] = '\0';
			}

			else {
				if (seq[0] == '\0' && end == 0) {		// blank line�� ���, �� seq �����̹Ƿ� motif ���� ��� ���
					whole[len] = '\0';		// seq ���� ���� (��ħǥ �����)	

					fprintf(fspace, "%d\n", len);
					fprintf(fspace, "%s\n", whole);

					fprintf(fsum2, "\t>Whole sequence length: %d\n", len);

					for (p = 0; p < LENGTH; p++)	// motif counter �迭 �ʱ�ȭ
						for (q = 0; q < AATYPE; q++)
							for (r = 0; r < AATYPE; r++)
								motif_db[p][q][r] = 0;

					motif_type = 0;					// motif ���� �ʱ�ȭ

					for (p = 0; p < LENGTH; p++) {	// motif counter ���
						for (q = 0; q < AATYPE; q++) {
							for (r = 0; r < AATYPE; r++) {
								if (p == 0) {
									cmp1[0] = '*';
									cmp1[1] = 'A' + q;
									cmp1[2] = 'A' + r;
								}

								if (p == 1) {
									cmp1[0] = 'A' + q;
									cmp1[1] = '*';
									cmp1[2] = 'A' + r;
								}

								if (p == 2) {
									cmp1[0] = 'A' + q;
									cmp1[1] = 'A' + r;
									cmp1[2] = '*';
								}

								for (m = 0; m <= len - LENGTH; m++) {
									for (n = 0; n < LENGTH; n++) {
										if (cmp1[n] == '*')
											cmp2[n] = '*';
										else
											cmp2[n] = whole[m + n];
									}

									if (strcmp(cmp1, cmp2) == 0) {
										if (motif_db[p][q][r] == 0)	// MOTIF ������ ���
											motif_type++;	

										motif_db[p][q][r]++;
									}
								}
							}
						}
					}

					fprintf(fout, "%d\n", motif_type);		// ��� ������ ù ���� whole sequence length, motif ����
					fprintf(fsum1, "\t>Types of trimer motifs: %d\n", motif_type);

					for (p = 0; p < LENGTH; p++) {		// motif counter ��� ��� ���
						for (q = 0; q < AATYPE; q++) {
							for (r = 0; r < AATYPE; r++) {
								if (motif_db[p][q][r] > 0) {	// whole seq ���� �����ϴ� motif�� ���!
									if (p == 0)
										fprintf(fout, "*%c%c %d\n", q + 'A', r + 'A', motif_db[p][q][r]);
									if (p == 1)
										fprintf(fout, "%c*%c %d\n", q + 'A', r + 'A', motif_db[p][q][r]);
									if (p == 2)
										fprintf(fout, "%c%c* %d\n", q + 'A', r + 'A', motif_db[p][q][r]);
								}
							}
						}
					}

					end = 1;				// blank line 2�� ���� �߻� �� ���� ó�� ����
				}

				else {						// general line�� ���, ���� �Է� �ޱ�
					for (m = 0; seq[m] != '\0'; m++) {
						whole[len] = seq[m];
						len++;
					}
				}
			}
		}

		fclose(fseq);
		fclose(fout);
		fclose(fspace);
	}

	fclose(ffile1);
	fclose(ffile2);
	fclose(fsum1);
	fclose(fsum2);
	fclose(fin);
}