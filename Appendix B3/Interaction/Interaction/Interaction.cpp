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
	FILE *fmotif_list, *fspace_list, *fmotif, *fspace, *fint, *fsum, *fname;
	char filename[NAME], seqname[NAME], newdir[NAME], whole[MAX], **motif, interaction[LENGTH * 2 + 1];
	int i, j, k, l, m, temp, len, motif_type, cnt, flag, total_int, total_type, num_proteins;

	fmotif_list = fopen("Raw_motifs/filename.txt", "r");
	fspace_list = fopen("Without_space/filename.txt", "r");

	system("mkdir Interaction_DB\\");
	fsum = fopen("Interaction_DB/summary.txt", "w");
	fname = fopen("Interaction_DB/filename.txt", "w");
	
	fprintf(fsum, "Summary of Interaction DB Construction:\n");

	while (1) {
		// Motif 파일 리스트 가져오기
		for (i = 0; i < NAME; i++) {
			filename[i] = '\0';
			newdir[i] = '\0';
		}

		fgets(filename, NAME, fmotif_list);
		fscanf(fmotif_list, "%d\n", &num_proteins);

		for (i = 0; i < NAME; i++) {
			if (filename[i] == 10) {
				filename[i] = '\0';
				break;
			}
		}

		if (i == 200)	break;		// 공백 라인이 생기면 가져오기 종료!

		strcat_s(newdir, NAME, "Raw_motifs/");
		strcat_s(newdir, NAME, filename);
		fmotif = fopen(newdir, "r");

		// 공백 없는 Seq 파일 리스트 가져오기
		for (i = 0; i < NAME; i++) {
			filename[i] = '\0';
			newdir[i] = '\0';
		}
		
		fgets(filename, NAME, fspace_list);
		fscanf(fspace_list, "%d\n", &num_proteins);
		filename[strlen(filename) - 1] = '\0';
		strcat_s(newdir, NAME, "Without_space/");
		strcat_s(newdir, NAME, filename);
		fspace = fopen(newdir, "r");
		
		// Interaction DB 파일 쓰기 시작
		for (i = 0; i < NAME; i++)
			newdir[i] = '\0';
		strcat_s(newdir, NAME, "Interaction_DB/");
		filename[strlen(filename) - 5] = 'i';
		filename[strlen(filename) - 4] = 'n';
		filename[strlen(filename) - 3] = 't';
		filename[strlen(filename) - 2] = 'D';
		filename[strlen(filename) - 1] = 'B';
		fprintf(fname, "%s\n", filename);	// Interaction DB 파일 목록 업데이트!
		strcat_s(newdir, NAME, filename);
		fint = fopen(newdir, "w");
		
		fprintf(fsum, "\n");
		fprintf(fname, "%d\n", num_proteins);

		while (1) {
			total_int = 0;
			total_type = 0;

			// Seqname 입력 받기
			for (j = 0; j < NAME; j++)
				seqname[j] = '\0';
			fgets(seqname, NAME, fmotif);
			if (seqname[0] == '\0')	break;	// 더이상 seqname 없으면 종료!

			for (j = 0; j < NAME; j++)
				seqname[j] = '\0';
			fgets(seqname, NAME, fspace);
			fprintf(fsum, "%s", &seqname[1]);
			printf("%s", &seqname);

			// Seqname을 Interation_DB 파일 내에 출력
			fprintf(fint, "%s", seqname);

			// 공백 없는 서열: 길이와 Whole Sequence 입력 받기
			fscanf(fspace, "%d\n", &len);
			for (j = 0; j < MAX; j++)
				whole[j] = '\0';
			fscanf(fspace, "%s\n\n", whole);

			// Motif 서열 : Motif 종류 갯수와 각 종류의 서열 입력 받기
			fscanf(fmotif, "%d\n", &motif_type);
			motif = (char **)malloc(sizeof(char*) * motif_type);
			for (j = 0; j < motif_type; j++)
				motif[j] = (char *)malloc(sizeof(char) * (LENGTH + 1));
			for (j = 0; j < motif_type; j++) {
				for (k = 0; k <= LENGTH; k++)
					motif[j][k] = '\0';
			}
			for (j = 0; j < motif_type; j++)
				fscanf(fmotif, "%s %d\n", motif[j], &temp);		// 서열만 입력받고 갯수는 temp로 처리해서 버림

			// Interaction 횟수 계산해서 출력
			for (j = 0; j < motif_type; j++) {
				for (k = 0; k < motif_type; k++) {
					for (l = 0 ; l < LENGTH; l++)
						interaction[l] = motif[j][l];
					for (l = LENGTH; l < LENGTH * 2; l++)
						interaction[l] = motif[k][l - LENGTH];
					interaction[LENGTH * 2] = '\0';

					cnt = 0;
					for (l = 0; l <= len - (LENGTH * 2); l++) {
						flag = 0;	// 0일 경우 같음, 1일 경우 다름
						for (m = 0; m < LENGTH * 2; m++) {
							if ((interaction[m] != whole[l + m]) && (interaction[m] != '*')) {
								flag = 1;
								break;
							}
						}

						if (flag == 0)	// 같다면 카운트
							cnt++;
					}

					if (cnt > 0) {
						fprintf(fint, "%s %d\n", interaction, cnt);
						total_int += cnt;
						total_type++;
					}
				}
			}

			fprintf(fsum, "\t>Total %d interaction(s) of %d interaction type(s)!\n", total_int, total_type);
			fprintf(fname, "%d ", total_type);	// Interaction DB 파일 목록 업데이트!

			for (j = 0; j < motif_type; j++)
				free(motif[j]);
			free(motif);
		}

		fprintf(fname, "\n");
		fclose(fmotif);
		fclose(fspace);
		fclose(fint);
	}

	fclose(fmotif_list);
	fclose(fspace_list);

	fclose(fsum);
	fclose(fname);
}